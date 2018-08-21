#include "spatial_discretization/finite_element_method/04_rhs.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"


namespace SpatialDiscretization
{

// 1D rhs
template<typename QuadratureType, typename Term>
void FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<1>, Equation::hasLaplaceOperator<Term>>::
transferRhsToWeakForm()
{
  LOG(TRACE)<<"transferRhsToWeakForm (1D)";

  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // get settings values
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  element_no_t nElements = mesh->nElementsLocal();
  double elementLength = mesh->meshWidth();
  dof_no_t nUnknownsLocal = mesh->nNodesLocalWithGhosts();

  LOG(DEBUG) << "Use settings nElements="<<nElements<<", elementLength="<<elementLength;

  // multiply factor to rhs
  // rhs *= stencil * elementLength
  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();

  // merge local changes on the vector
  rightHandSide.startVectorManipulation();
  
  // stencil values
  // stencil in 1D: 1/6*[1 _4_ 1] (element contribution: 1/6*[_2_ 1])
  const int center = 1;
  const double stencilCenter[3] = {1./6., 4./6., 1./6.};
  const double stencilSide[2] = {2./6., 1./6.};

  // get all entries
  std::vector<double> vectorValues;
  rightHandSide.getValuesWithGhosts(vectorValues, false);

  rightHandSide.zeroEntries();
  
  // loop over all dofs and set values with stencilCenter
  for (node_no_t dofNo = 1; dofNo < nUnknownsLocal-1; dofNo++)
  {
    double value =
      (stencilCenter[center-1]*vectorValues[dofNo-1]
      + stencilCenter[center]*vectorValues[dofNo]
      + stencilCenter[center+1]*vectorValues[dofNo+1]) * elementLength;

    rightHandSide.setValue(dofNo, value, INSERT_VALUES);
  }
  
  rightHandSide.finishVectorManipulation();
  rightHandSide.startVectorManipulation();
  rightHandSide.getValuesWithGhosts(vectorValues);
  
  // set values for boundaries with stencilSide
  node_no_t dofNo = 0;
  double value =
    (stencilSide[0]*vectorValues[dofNo]
    + stencilSide[1]*vectorValues[dofNo+1]) * elementLength;
  rightHandSide.setValue(0, value, ADD_VALUES);

  dofNo = nUnknownsLocal-1;
  value =
    (stencilSide[0]*vectorValues[dofNo]
    + stencilSide[1]*vectorValues[dofNo-1]) * elementLength;
  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // call VecAssemblyBegin, VecAssemblyEnd
  rightHandSide.finishVectorManipulation();
}

// 2D rhs
template<typename QuadratureType, typename Term>
void FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<2>, Equation::hasLaplaceOperator<Term>>::
transferRhsToWeakForm()
{
  LOG(TRACE)<<"transferRhsToWeakForm (2D)";

  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // get settings values
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  element_no_t nElements0 = mesh->nElementsPerCoordinateDirectionLocal(0);
  element_no_t nElements1 = mesh->nElementsPerCoordinateDirectionLocal(1);
  node_no_t nNodes0 = mesh->nNodesLocalWithGhosts(0);
  node_no_t nNodes1 = mesh->nNodesLocalWithGhosts(1);
  double elementLength0 = mesh->meshWidth();
  double elementLength1 = mesh->meshWidth();
  double integralFactor = elementLength0*elementLength1;

  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();

  // merge local changes on the vector
  rightHandSide.startVectorManipulation();
  
  // stencil values

  // stencil for rhs in 2D:      [1  4   1] (element contribution:      [ 2  1])
  //                        1/36*[4 _16_ 4]                        1/36*[_4_ 2]
  //                             [1  4   1]

  const int center = 1;
  const double stencilCenter[3][3] = {
    {1./36, 4./36,  1./36},
    {4./36, 16./36, 4./36},
    {1./36, 4./36,  1./36}};

  const double stencilEdge[2][3] = {
    {1./36, 4./36, 1./36},
    {2./36, 8./36, 2./36}
  };

  const double stencilCorner[2][2] = {
    {1./36, 2./36},
    {2./36, 4./36}
  };

  auto dofIndex = [&nNodes0](int x, int y){return y*nNodes0 + x;};
  double value;

  // get all values
  std::vector<double> vectorValues;
  rightHandSide.getValuesWithGhosts(vectorValues);

  rightHandSide.zeroEntries();
  // loop over all dofs and set values with stencilCenter
  // set entries for interior nodes
  for (int y=1; y<nNodes1-1; y++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      value = 0;
      for (int i=-1; i<=1; i++) // x
      {
        for (int j=-1; j<=1; j++) // y
        {
          value += stencilCenter[center+i][center+j] * vectorValues[dofIndex(x+i, y+j)];
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofIndex(x,y), value, INSERT_VALUES);
    }
  }

  rightHandSide.finishVectorManipulation();
  rightHandSide.startVectorManipulation();
  rightHandSide.getValuesWithGhosts(vectorValues);
  
  // set entries for boundary nodes on edges
  // left boundary (x=0)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    node_no_t dofNo = dofIndex(x,y);

    value = 0;
    for (int i=-1; i<=0; i++) // -x
    {
      for (int j=-1; j<=1; j++) // y
      {
        value += stencilEdge[center+i][center+j] * vectorValues[dofIndex(x-i,y+j)];
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // right boundary (x=nNodes0-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    node_no_t dofNo = dofIndex(x,y);

    value = 0;
    for (int i=-1; i<=0; i++) // x
    {
      for (int j=-1; j<=1; j++) // y
      {
        value += stencilEdge[center+i][center+j] * vectorValues[dofIndex(x+i,y+j)];
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // bottom boundary (y=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    node_no_t dofNo = dofIndex(x,y);

    value = 0;
    for (int i=-1; i<=1; i++) // x
    {
      for (int j=-1; j<=0; j++) // -y
      {
        value += stencilEdge[center+j][center+i] * vectorValues[dofIndex(x+i,y-j)];
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // top boundary (y=nNodes1-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    node_no_t dofNo = dofIndex(x,y);

    value = 0;
    for (int i=-1; i<=1; i++) // x
    {
      for (int j=-1; j<=0; j++) // y
      {
        value += stencilEdge[center+j][center+i] * vectorValues[dofIndex(x+i,y+j)];
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // corner nodes
  int x,y;
  node_no_t dofNo;

  // bottom left (x=0, y=0)
  x = 0;
  y = 0;
  dofNo = dofIndex(x,y);

  value = 0;
  for (int i=-1; i<=0; i++) // -x
  {
    for (int j=-1; j<=0; j++) // -y
    {
      value += stencilCorner[center+i][center+j] * vectorValues[dofIndex(x-i,y-j)];
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // bottom right (x=nNodes0-1, y=0)
  x = nNodes0-1;
  y = 0;
  dofNo = dofIndex(x,y);

  value = 0;
  for (int i=-1; i<=0; i++) // x
  {
    for (int j=-1; j<=0; j++) // -y
    {
      value += stencilCorner[center+i][center+j] * vectorValues[dofIndex(x+i,y-j)];
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // top left (x=0, y=nNodes1-1)
  x = 0;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);

  value = 0;
  for (int i=-1; i<=0; i++) // -x
  {
    for (int j=-1; j<=0; j++) // y
    {
      value += stencilCorner[center+i][center+j] * vectorValues[dofIndex(x-i,y+j)];
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // top right (x=nNodes0-1, y=nNodes1-1)
  x = nNodes0-1;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);

  value = 0;
  for (int i=-1; i<=0; i++) // x
  {
    for (int j=-1; j<=0; j++) // y
    {
      value += stencilCorner[center+i][center+j] * vectorValues[dofIndex(x+i,y+j)];
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // call VecAssemblyBegin, VecAssemblyEnd
  rightHandSide.finishVectorManipulation();
}

// 3D rhs
template<typename QuadratureType, typename Term>
void FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<3>, Equation::hasLaplaceOperator<Term>>::
transferRhsToWeakForm()
{
  LOG(TRACE)<<"transferRhsToWeakForm (3D)";

  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // get settings values
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  element_no_t nElements0 = mesh->nElementsPerCoordinateDirectionLocal(0);
  element_no_t nElements1 = mesh->nElementsPerCoordinateDirectionLocal(1);
  element_no_t nElements2 = mesh->nElementsPerCoordinateDirectionLocal(2);
  node_no_t nNodes0 = mesh->nNodesLocalWithGhosts(0);
  node_no_t nNodes1 = mesh->nNodesLocalWithGhosts(1);
  node_no_t nNodes2 = mesh->nNodesLocalWithGhosts(2);
  double elementLength0 = mesh->meshWidth();
  double elementLength1 = mesh->meshWidth();
  double elementLength2 = mesh->meshWidth();
  double integralFactor = elementLength0*elementLength1*elementLength2;

  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();

  // merge local changes on the vector
  rightHandSide.startVectorManipulation();
  
  // stencil values

  // stencil for rhs in 3D: bottom: [1  4   1]   (element contribution:  center:[ 4  2]
  //                          1/216*[4 _16_ 4]                            1/216*[_8_ 4]
  //                                [1  4   1]
  //                                                                     bottom:[ 2  1]
  //                        center: [ 4  16   4]                          1/216*[ 4  2] )
  //                          1/216*[16 _64_ 16]
  //                                [ 4  16   4]
  //
  //                        top: like bottom

  // coordinate system
  // x axis: left -> right
  // y axis: front -> back
  // z axis: bottom -> top


  const int center = 1;
  const double stencilCenter[3][3][3] = {
    {{1./216, 4./216,  1./216},   //bottom
    {4./216, 16./216, 4./216},
    {1./216, 4./216,  1./216}},
    {{4./216, 16./216,  4./216},   //center
    {16./216, 64./216, 16./216},
    {4./216, 16./216,  4./216}},
    {{1./216, 4./216,  1./216},   //top
    {4./216, 16./216, 4./216},
    {1./216, 4./216,  1./216}}
  };

  const double stencilBoundarySurface[2][3][3] = {
    {{1./216, 4./216,  1./216},   //bottom
    {4./216, 16./216, 4./216},
    {1./216, 4./216,  1./216}},
    {{2./216, 8./216,  2./216},   //center
    {8./216, 32./216, 8./216},
    {2./216, 8./216,  2./216}}
  };
  const double stencilBoundaryEdge[2][2][3] = {
    {{1./216, 4./216, 1./216},
    {2./216, 8./216, 2./216}},    //bottom
    {{2./216, 8./216, 2./216},
    {4./216, 16./216, 4./216}}    //center
  };

  const double stencilCorner[2][2][2] = {
    {{1./216, 2./216},
    {2./216, 4./216}},    //bottom
    {{2./216, 4./216},
    {4./216, 8./216}},    //center
  };

  auto dofIndex = [&nNodes0, &nNodes1](int x, int y, int z){
    return z*nNodes0*nNodes1 + y*nNodes0 + x;
  };
  double value;

  // get all values
  std::vector<double> vectorValues;
  rightHandSide.getValuesWithGhosts(vectorValues);

  rightHandSide.zeroEntries();
  // loop over all dofs and set values with stencilCenter
  // set entries for interior nodes
  for (int z=1; z<nNodes2-1; z++)
  {
    for (int y=1; y<nNodes1-1; y++)
    {
      for (int x=1; x<nNodes0-1; x++)
      {
        value = 0;
        for (int i=-1; i<=1; i++) // x
        {
          for (int j=-1; j<=1; j++) // y
          {
            for (int k=-1; k<=1; k++) // z
            {
              value += stencilCenter[center+i][center+j][center+k] * vectorValues[dofIndex(x+i, y+j, z+k)];
            }
          }
        }
        value *= integralFactor;

        rightHandSide.setValue(dofIndex(x,y,z), value, INSERT_VALUES);
      }
    }
  }

  rightHandSide.finishVectorManipulation();
  rightHandSide.startVectorManipulation();
  rightHandSide.getValuesWithGhosts(vectorValues);
  
  // set entries for boundary nodes on surface boundaries
  // left boundary (x = 0)
  for (int z=1; z<nNodes2-1; z++)
  {
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = 0;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=0; i++)    // -x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value += stencilBoundarySurface[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y+j,z+k)];
          }
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofNo, value, ADD_VALUES);
    }
  }

  // right boundary (x = nNodes0-1)
  for (int z=1; z<nNodes2-1; z++)
  {
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = nNodes0-1;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=0; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value += stencilBoundarySurface[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y+j,z+k)];
          }
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofNo, value, ADD_VALUES);
    }
  }

  // front boundary (y = 0)
  for (int z=1; z<nNodes2-1; z++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = 0;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // -y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value += stencilBoundarySurface[center+j][center+i][center+k] * vectorValues[dofIndex(x+i,y-j,z+k)];
          }
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofNo, value, ADD_VALUES);
    }
  }

  // back boundary (y = nNodes1-1)
  for (int z=1; z<nNodes2-1; z++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = nNodes1-1;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value += stencilBoundarySurface[center+j][center+i][center+k] * vectorValues[dofIndex(x+i,y+j,z+k)];
          }
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofNo, value, ADD_VALUES);
    }
  }

  // bottom boundary (z = 0)
  for (int y=1; y<nNodes1-1; y++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      int z = 0;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // -z
          {
            value += stencilBoundarySurface[center+k][center+i][center+j] * vectorValues[dofIndex(x+i,y+j,z-k)];
          }
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofNo, value, ADD_VALUES);
    }
  }

  // top boundary (z = nNodes2-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      int z = nNodes2-1;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // z
          {
            value += stencilBoundarySurface[center+k][center+i][center+j] * vectorValues[dofIndex(x+i,y+j,z+k)];
          }
        }
      }
      value *= integralFactor;

      rightHandSide.setValue(dofNo, value, ADD_VALUES);
    }
  }

  // set entries for boundary nodes on edge boundaries
  // bottom left (x=0,z=0)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    int z = 0;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value += stencilBoundaryEdge[center+i][center+k][center+j] * vectorValues[dofIndex(x-i,y+j,z-k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // bottom right (x=nNodes0-1,z=0)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    int z = 0;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value += stencilBoundaryEdge[center+i][center+k][center+j] * vectorValues[dofIndex(x+i,y+j,z-k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // top left (x=0,z=nNodes2-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    int z = nNodes2-1;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value += stencilBoundaryEdge[center+i][center+k][center+j] * vectorValues[dofIndex(x-i,y+j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // top right (x=nNodes0-1,z=nNodes2-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    int z = nNodes2-1;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value += stencilBoundaryEdge[center+i][center+k][center+j] * vectorValues[dofIndex(x+i,y+j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // bottom front (y=0,z=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    int z = 0;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value += stencilBoundaryEdge[center+j][center+k][center+i] * vectorValues[dofIndex(x+i,y-j,z-k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // bottom back (y=nNodes1-1,z=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    int z = 0;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value += stencilBoundaryEdge[center+j][center+k][center+i] * vectorValues[dofIndex(x+i,y+j,z-k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // top front (y=0,z=nNodes2-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    int z = nNodes2-1;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value += stencilBoundaryEdge[center+j][center+k][center+i] * vectorValues[dofIndex(x+i,y-j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // top back (y=nNodes1-1,z=nNodes2-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    int z = nNodes2-1;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value += stencilBoundaryEdge[center+j][center+k][center+i] * vectorValues[dofIndex(x+i,y+j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // left front (x=0,y=0)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = 0;
    int y = 0;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value += stencilBoundaryEdge[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y-j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // left back (x=0,y=nNodes1-1)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = 0;
    int y = nNodes1-1;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value += stencilBoundaryEdge[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y+j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // right front (x=nNodes0-1,y=0)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = nNodes0-1;
    int y = 0;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value += stencilBoundaryEdge[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y-j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // right back (x=nNodes0-1,y=nNodes1-1)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = nNodes0-1;
    int y = nNodes1-1;
    node_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value += stencilBoundaryEdge[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y+j,z+k)];
        }
      }
    }
    value *= integralFactor;

    rightHandSide.setValue(dofNo, value, ADD_VALUES);
  }

  // corner nodes
  int x,y,z;
  node_no_t dofNo;

  // bottom front left (x=0,y=0,z=0)
  x = 0;
  y = 0;
  z = 0;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // -x
  {
    for (int j=-1; j<=0; j++)   // -y
    {
      for (int k=-1; k<=0; k++)   // -z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y-j,z-k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // bottom front right (x=nNodes0-1,y=0,z=0)
  x = nNodes0-1;
  y = 0;
  z = 0;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // x
  {
    for (int j=-1; j<=0; j++)   // -y
    {
      for (int k=-1; k<=0; k++)   // -z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y-j,z-k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // bottom back left (x=0,y=nNodes1-1,z=0)
  x = 0;
  y = nNodes1-1;
  z = 0;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // -x
  {
    for (int j=-1; j<=0; j++)   // y
    {
      for (int k=-1; k<=0; k++)   // -z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y+j,z-k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // bottom back right (x=nNodes0-1,y=nNodes1-1,z=0)
  x = nNodes0-1;
  y = nNodes1-1;
  z = 0;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // x
  {
    for (int j=-1; j<=0; j++)   // y
    {
      for (int k=-1; k<=0; k++)   // -z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y+j,z-k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // top front left (x=0,y=0,z=nNodes2-1)
  x = 0;
  y = 0;
  z = nNodes2-1;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // -x
  {
    for (int j=-1; j<=0; j++)   // -y
    {
      for (int k=-1; k<=0; k++)   // z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y-j,z+k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // top front right (x=nNodes0-1,y=0,z=nNodes2-1)
  x = nNodes0-1;
  y = 0;
  z = nNodes2-1;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // x
  {
    for (int j=-1; j<=0; j++)   // -y
    {
      for (int k=-1; k<=0; k++)   // z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y-j,z+k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // top back left (x=0,y=nNodes1-1,z=nNodes2-1)
  x = 0;
  y = nNodes1-1;
  z = nNodes2-1;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // -x
  {
    for (int j=-1; j<=0; j++)   // y
    {
      for (int k=-1; k<=0; k++)   // z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x-i,y+j,z+k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // top back right (x=nNodes0-1,y=nNodes1-1,z=nNodes2-1)
  x = nNodes0-1;
  y = nNodes1-1;
  z = nNodes2-1;
  dofNo = dofIndex(x,y,z);

  value = 0;
  for (int i=-1; i<=0; i++)    // x
  {
    for (int j=-1; j<=0; j++)   // y
    {
      for (int k=-1; k<=0; k++)   // z
      {
        value += stencilCorner[center+i][center+j][center+k] * vectorValues[dofIndex(x+i,y+j,z+k)];
      }
    }
  }
  value *= integralFactor;

  rightHandSide.setValue(dofNo, value, ADD_VALUES);

  // call VecAssemblyBegin, VecAssemblyEnd
  rightHandSide.finishVectorManipulation();
}

// 1D massMatrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<1>, Equation::hasLaplaceOperator<Term>>::
setMassMatrix()
{
  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // check if matrix discretization matrix exists
  if (!this->data_.massMatrixInitialized())
  {
    this->data_.initializeMassMatrix();

    // set entries of matrix
    LOG(DEBUG)<<"createMassMatrix 1D";

    // get settings values
    std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
    element_no_t nElements = mesh->nElementsLocal();
    double elementLength = mesh->meshWidth();

    dof_no_t nUnknownsLocal = mesh->nNodesLocalWithGhosts();

    LOG(DEBUG) << "Use settings nElements="<<nElements<<", elementLength="<<elementLength;

    // multiply factor to rhs
    // rhs *= stencil * elementLength
    std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix = this->data_.massMatrix();
    this->data_.massMatrix().zeroEntries();
    
    // massMatrix * f_strong = rhs_weak
    // row of massMatrix: contributions to a single entry in rhs_weak

    // stencil values
    // stencil in 1D: 1/6*[1 _4_ 1] (element contribution: 1/6*[_2_ 1])
    const int center = 1;
    const double stencilCenter[3] = {1./6.*1, 1./6.*4, 1./6.*1};
    const double stencilSide[2] = {1./6.*2, 1./6.*1};

    // loop over all dofs and set values in massMatrix with stencilCenter
    for (node_no_t dofNo = 1; dofNo < nUnknownsLocal-1; dofNo++)
    {
      massMatrix->setValue(dofNo, dofNo-1, stencilCenter[center-1] * elementLength, INSERT_VALUES);
      massMatrix->setValue(dofNo, dofNo,   stencilCenter[center]   * elementLength, INSERT_VALUES);
      massMatrix->setValue(dofNo, dofNo+1, stencilCenter[center+1] * elementLength, INSERT_VALUES);
    }
    
    // call MatAssemblyBegin, MatAssemblyEnd
    massMatrix->assembly(MAT_FLUSH_ASSEMBLY);

    // set values for boundaries with stencilSide
    node_no_t dofNo = 0;
    massMatrix->setValue(dofNo, dofNo,   stencilSide[0] * elementLength, ADD_VALUES);
    massMatrix->setValue(dofNo, dofNo+1, stencilSide[1] * elementLength, ADD_VALUES);

    dofNo = nUnknownsLocal-1;
    massMatrix->setValue(dofNo, dofNo,   stencilSide[0] * elementLength, ADD_VALUES);
    massMatrix->setValue(dofNo, dofNo-1, stencilSide[1] * elementLength, ADD_VALUES);

    // call MatAssemblyBegin, MatAssemblyEnd
    massMatrix->assembly(MAT_FINAL_ASSEMBLY);
  }
}

// 2D massMatrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<2>, Equation::hasLaplaceOperator<Term>>::
setMassMatrix()
{
  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // check if matrix discretization matrix exists
  if (!this->data_.massMatrixInitialized())
  {
    this->data_.initializeMassMatrix();

    // set entries of matrix
    LOG(DEBUG)<<"createMassMatrix 2D";

    // get settings values
    std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
    element_no_t nElements0 = mesh->nElementsPerCoordinateDirectionLocal(0);
    element_no_t nElements1 = mesh->nElementsPerCoordinateDirectionLocal(1);
    node_no_t nNodes0 = mesh->nNodesLocalWithGhosts(0);
    node_no_t nNodes1 = mesh->nNodesLocalWithGhosts(1);
    double elementLength0 = mesh->meshWidth();
    double elementLength1 = mesh->meshWidth();
    double integralFactor = elementLength0*elementLength1;

    LOG(DEBUG) << "elementLength0: " << elementLength0;

    // multiply factor to rhs
    // rhs *= stencil * elementLength
    std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix = this->data_.massMatrix();
    this->data_.massMatrix().zeroEntries();
    
    // massMatrix * f_strong = rhs_weak
    // row of massMatrix: contributions to a single entry in rhs_weak

    // stencil values

    // stencil for rhs in 2D:      [1  4   1] (element contribution:      [ 2  1])
    //                        1/36*[4 _16_ 4]                        1/36*[_4_ 2]
    //                             [1  4   1]

    const int center = 1;
    const double stencilCenter[3][3] = {
      {1./36, 4./36,  1./36},
      {4./36, 16./36, 4./36},
      {1./36, 4./36,  1./36}};

    const double stencilEdge[2][3] = {
      {2./36, 8./36, 2./36},
      {1./36, 4./36, 1./36}
    };

    const double stencilCorner[2][2] = {
      {4./36, 2./36},
      {2./36, 1./36}
    };

    auto dofIndex = [&nElements0, &nElements1](int x, int y){return y*(nElements0+1) + x;};

    // loop over all dofs and set values with stencilCenter
    // set entries for interior nodes
    for (int y=1; y<nNodes1-1; y++)
    {
      for (int x=1; x<nNodes0-1; x++)
      {
        node_no_t dofNo = dofIndex(x,y);
        for (int i=-1; i<=1; i++)
        {
          for (int j=-1; j<=1; j++)
          {
            node_no_t secondDofNo = dofIndex(x+j, y+i);
            double factor = stencilCenter[center+i][center+j] * integralFactor;
            massMatrix->setValue(dofNo, secondDofNo, factor, INSERT_VALUES);
          }
        }
      }
    }

    // call MatAssemblyBegin, MatAssemblyEnd
    massMatrix->assembly(MAT_FLUSH_ASSEMBLY);

    // set entries for boundary nodes on edges
    // left boundary
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = 0;
      node_no_t dofNo = dofIndex(x,y);

      for (int i=-1; i<=1; i++)
      {
        node_no_t secondDofNo = dofIndex(x, y+i);
        double factor = stencilEdge[0][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);

        secondDofNo = dofIndex(x+1, y+i);
        factor = stencilEdge[1][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);
      }
    }

    // right boundary
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = nNodes0-1;
      node_no_t dofNo = dofIndex(x,y);

      for (int i=-1; i<=1; i++)
      {
        node_no_t secondDofNo = dofIndex(x, y+i);
        double factor = stencilEdge[0][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);

        secondDofNo = dofIndex(x-1, y+i);
        factor = stencilEdge[1][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);
      }
    }

    // bottom boundary
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = 0;
      node_no_t dofNo = dofIndex(x,y);

      for (int i=-1; i<=1; i++)
      {
        node_no_t secondDofNo = dofIndex(x+i, y);
        double factor = stencilEdge[0][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);

        secondDofNo = dofIndex(x+i, y+1);
        factor = stencilEdge[1][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);
      }
    }

    // top boundary
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = nNodes1-1;
      node_no_t dofNo = dofIndex(x,y);

      for (int i=-1; i<=1; i++)
      {
        node_no_t secondDofNo = dofIndex(x+i, y);
        double factor = stencilEdge[0][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);

        secondDofNo = dofIndex(x+i, y-1);
        factor = stencilEdge[1][center+i] * integralFactor;
        massMatrix->setValue(dofNo, secondDofNo, factor, ADD_VALUES);
      }
    }

    // corner nodes
    // bottom left
    int x = 0;
    int y = 0;
    node_no_t dofNo = dofIndex(x,y);

    massMatrix->setValue(dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x+1,y),   stencilCorner[0][1] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x,y+1),   stencilCorner[1][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x+1,y+1), stencilCorner[1][1] * integralFactor, ADD_VALUES);

    // bottom right
    x = nNodes0-1;
    y = 0;
    dofNo = dofIndex(x,y);

    massMatrix->setValue(dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x-1,y),   stencilCorner[0][1] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x,y+1),   stencilCorner[1][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x-1,y+1), stencilCorner[1][1] * integralFactor, ADD_VALUES);

    // top left
    x = 0;
    y = nNodes1-1;
    dofNo = dofIndex(x,y);

    massMatrix->setValue(dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x+1,y),   stencilCorner[0][1] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x,y-1),   stencilCorner[1][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x+1,y-1), stencilCorner[1][1] * integralFactor, ADD_VALUES);

    // top right
    x = nNodes0-1;
    y = nNodes1-1;
    dofNo = dofIndex(x,y);

    massMatrix->setValue(dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x-1,y),   stencilCorner[0][1] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x,y-1),   stencilCorner[1][0] * integralFactor, ADD_VALUES);
    massMatrix->setValue(dofNo, dofIndex(x-1,y-1), stencilCorner[1][1] * integralFactor, ADD_VALUES);

    // call MatAssemblyBegin, MatAssemblyEnd
    massMatrix->assembly(MAT_FINAL_ASSEMBLY);
  }
}

// 3D massMatrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<3>, Equation::hasLaplaceOperator<Term>>::
setMassMatrix()
{
  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // check if matrix discretization matrix exists
  if (!this->data_.massMatrixInitialized())
  {
    this->data_.initializeMassMatrix();

    // set entries of matrix
    LOG(DEBUG)<<"createMassMatrix 3D";

    // get settings values
    std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
    element_no_t nElements0 = mesh->nElementsPerCoordinateDirectionLocal(0);
    element_no_t nElements1 = mesh->nElementsPerCoordinateDirectionLocal(1);
    element_no_t nElements2 = mesh->nElementsPerCoordinateDirectionLocal(2);
    node_no_t nNodes0 = mesh->nNodesLocalWithGhosts(0);
    node_no_t nNodes1 = mesh->nNodesLocalWithGhosts(1); 
    node_no_t nNodes2 = mesh->nNodesLocalWithGhosts(2);
    double elementLength0 = mesh->meshWidth();
    double elementLength1 = mesh->meshWidth();
    double elementLength2 = mesh->meshWidth();
    double integralFactor = elementLength0*elementLength1*elementLength2;

    // multiply factor to rhs
    // rhs *= stencil * elementLength
    std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix = this->data_.massMatrix();
    this->data_.massMatrix().zeroEntries();
    
    // massMatrix * f_strong = rhs_weak
    // row of massMatrix: contributions to a single entry in rhs_weak

    // stencil values

    // stencil for rhs in 3D: bottom: [1  4   1]   (element contribution:  center:[ 4  2]
    //                          1/216*[4 _16_ 4]                            1/216*[_8_ 4]
    //                                [1  4   1]
    //                                                                     bottom:[ 2  1]
    //                        center: [ 4  16   4]                          1/216*[ 4  2] )
    //                          1/216*[16 _64_ 16]
    //                                [ 4  16   4]
    //
    //                        top: like bottom

    // coordinate system
    // x axis: left -> right
    // y axis: front -> back
    // z axis: bottom -> top


    const int center = 1;
    const double stencilCenter[3][3][3] = {
      {{1./216, 4./216,  1./216},   //bottom
      {4./216, 16./216, 4./216},
      {1./216, 4./216,  1./216}},
      {{4./216, 16./216,  4./216},   //center
      {16./216, 64./216, 16./216},
      {4./216, 16./216,  4./216}},
      {{1./216, 4./216,  1./216},   //top
      {4./216, 16./216, 4./216},
      {1./216, 4./216,  1./216}}
    };

    const double stencilBoundarySurface[2][3][3] = {
      {{1./216, 4./216,  1./216},   //bottom
      {4./216, 16./216, 4./216},
      {1./216, 4./216,  1./216}},
      {{2./216, 8./216,  2./216},   //center
      {8./216, 32./216, 8./216},
      {2./216, 8./216,  2./216}}
    };
    const double stencilBoundaryEdge[2][2][3] = {
      {{1./216, 4./216, 1./216},
      {2./216, 8./216, 2./216}},    //bottom
      {{2./216, 8./216, 2./216},
      {4./216, 16./216, 4./216}}    //center
    };

    const double stencilCorner[2][2][2] = {
      {{1./216, 2./216},
      {2./216, 4./216}},    //bottom
      {{2./216, 4./216},
      {4./216, 8./216}},    //center
    };

    auto dofIndex = [&nNodes0, &nNodes1](int x, int y, int z){
      return z*nNodes0*nNodes1 + y*nNodes0 + x;
    };

    double value;
    node_no_t dofNo;

    // loop over all dofs and set values with stencilCenter
    // set entries for interior nodes
    for (int z=1; z<nNodes2-1; z++)
    {
      for (int y=1; y<nNodes1-1; y++)
      {
        for (int x=1; x<nNodes0-1; x++)
        {
          dofNo = dofIndex(x, y, z);
          for (int i=-1; i<=1; i++) // x
          {
            for (int j=-1; j<=1; j++) // y
            {
              for (int k=-1; k<=1; k++) // z
              {
                value = stencilCenter[center+i][center+j][center+k]*integralFactor;
                //                 matrix           row    column
                massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES);
              }
            }
          }
        }
      }
    }

    // call MatAssemblyBegin, MatAssemblyEnd
    massMatrix->assembly(MAT_FLUSH_ASSEMBLY);

    // set entries for boundary nodes on surface boundaries
    // left boundary (x = 0)
    for (int z=1; z<nNodes2-1; z++)
    {
      for (int y=1; y<nNodes1-1; y++)
      {
        int x = 0;
        node_no_t dofNo = dofIndex(x,y,z);
        for (int i=-1; i<=0; i++)    // -x
        {
          for (int j=-1; j<=1; j++)   // y
          {
            for (int k=-1; k<=1; k++)   // z
            {
              value = stencilBoundarySurface[center+i][center+j][center+k]*integralFactor;
              //                 matrix           row    column
              massMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
            }
          }
        }
      }
    }

    // right boundary (x = nNodes0-1)
    for (int z=1; z<nNodes2-1; z++)
    {
      for (int y=1; y<nNodes1-1; y++)
      {
        int x = nNodes0-1;
        node_no_t dofNo = dofIndex(x,y,z);
        for (int i=-1; i<=0; i++)    // x
        {
          for (int j=-1; j<=1; j++)   // y
          {
            for (int k=-1; k<=1; k++)   // z
            {
              value = stencilBoundarySurface[center+i][center+j][center+k]*integralFactor;
              //                 matrix           row    column
              massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
            }
          }
        }
      }
    }

    // front boundary (y = 0)
    for (int z=1; z<nNodes2-1; z++)
    {
      for (int x=1; x<nNodes0-1; x++)
      {
        int y = 0;
        node_no_t dofNo = dofIndex(x,y,z);
        for (int i=-1; i<=1; i++)    // x
        {
          for (int j=-1; j<=0; j++)   // -y
          {
            for (int k=-1; k<=1; k++)   // z
            {
              value = stencilBoundarySurface[center+j][center+i][center+k]*integralFactor;
              //                 matrix           row    column
              massMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
            }
          }
        }
      }
    }

    // back boundary (y = nNodes1-1)
    for (int z=1; z<nNodes2-1; z++)
    {
      for (int x=1; x<nNodes0-1; x++)
      {
        int y = nNodes1-1;
        node_no_t dofNo = dofIndex(x,y,z);
        for (int i=-1; i<=1; i++)    // x
        {
          for (int j=-1; j<=0; j++)   // y
          {
            for (int k=-1; k<=1; k++)   // z
            {
              value = stencilBoundarySurface[center+j][center+i][center+k]*integralFactor;
              //                 matrix           row    column
              massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
            }
          }
        }
      }
    }

    // bottom boundary (z = 0)
    for (int y=1; y<nNodes1-1; y++)
    {
      for (int x=1; x<nNodes0-1; x++)
      {
        int z = 0;
        node_no_t dofNo = dofIndex(x,y,z);
        for (int i=-1; i<=1; i++)    // x
        {
          for (int j=-1; j<=1; j++)   // y
          {
            for (int k=-1; k<=0; k++)   // -z
            {
              value = stencilBoundarySurface[center+k][center+i][center+j]*integralFactor;
              //                 matrix           row    column
              massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
            }
          }
        }
      }
    }

    // top boundary (z = nNodes2-1)
    for (int y=1; y<nNodes1-1; y++)
    {
      for (int x=1; x<nNodes0-1; x++)
      {
        int z = nNodes2-1;
        node_no_t dofNo = dofIndex(x,y,z);
        for (int i=-1; i<=1; i++)    // x
        {
          for (int j=-1; j<=1; j++)   // y
          {
            for (int k=-1; k<=0; k++)   // z
            {
              value = stencilBoundarySurface[center+k][center+i][center+j]*integralFactor;
              //                 matrix           row    column
              massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
            }
          }
        }
      }
    }

    // set entries for boundary nodes on edge boundaries
    // bottom left (x=0,z=0)
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = 0;
      int z = 0;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // -x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // -z
          {
            value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x-i, y+j, z-k), value, ADD_VALUES);
          }
        }
      }
    }

    // bottom right (x=nNodes0-1,z=0)
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = nNodes0-1;
      int z = 0;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // -z
          {
            value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
          }
        }
      }
    }

    // top left (x=0,z=nNodes2-1)
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = 0;
      int z = nNodes2-1;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // -x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // z
          {
            value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // top right (x=nNodes0-1,z=nNodes2-1)
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = nNodes0-1;
      int z = nNodes2-1;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // z
          {
            value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // bottom front (y=0,z=0)
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = 0;
      int z = 0;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // -y
        {
          for (int k=-1; k<=0; k++)   // -z
          {
            value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y-j, z-k), value, ADD_VALUES);
          }
        }
      }
    }

    // bottom back (y=nNodes1-1,z=0)
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = nNodes1-1;
      int z = 0;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // -z
          {
            value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
          }
        }
      }
    }

    // top front (y=0,z=nNodes2-1)
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = 0;
      int z = nNodes2-1;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // -y
        {
          for (int k=-1; k<=0; k++)   // z
          {
            value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // top back (y=nNodes1-1,z=nNodes2-1)
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = nNodes1-1;
      int z = nNodes2-1;
      node_no_t dofNo = dofIndex(x,y,z);

      value = 0;
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // z
          {
            value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // left front (x=0,y=0)
    for (int z=1; z<nNodes2-1; z++)
    {
      int x = 0;
      int y = 0;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // -x
      {
        for (int j=-1; j<=0; j++)   // -y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x-i, y-j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // left back (x=0,y=nNodes1-1)
    for (int z=1; z<nNodes2-1; z++)
    {
      int x = 0;
      int y = nNodes1-1;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // -x
      {
        for (int j=-1; j<=0; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // right front (x=nNodes0-1,y=0)
    for (int z=1; z<nNodes2-1; z++)
    {
      int x = nNodes0-1;
      int y = 0;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // -y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // right back (x=nNodes0-1,y=nNodes1-1)
    for (int z=1; z<nNodes2-1; z++)
    {
      int x = nNodes0-1;
      int y = nNodes1-1;
      node_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
            //                 matrix           row    column
            massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
          }
        }
      }
    }

    // corner nodes
    int x,y,z;

    // bottom front left (x=0,y=0,z=0)
    x = 0;
    y = 0;
    z = 0;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x-i, y-j, z-k), value, ADD_VALUES);
        }
      }
    }

    // bottom front right (x=nNodes0-1,y=0,z=0)
    x = nNodes0-1;
    y = 0;
    z = 0;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x+i, y-j, z-k), value, ADD_VALUES);
        }
      }
    }

    // bottom back left (x=0,y=nNodes1-1,z=0)
    x = 0;
    y = nNodes1-1;
    z = 0;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x-i, y+j, z-k), value, ADD_VALUES);
        }
      }
    }

    // bottom back right (x=nNodes0-1,y=nNodes1-1,z=0)
    x = nNodes0-1;
    y = nNodes1-1;
    z = 0;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
        }
      }
    }

    // top front left (x=0,y=0,z=nNodes2-1)
    x = 0;
    y = 0;
    z = nNodes2-1;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x-i, y-j, z+k), value, ADD_VALUES);
        }
      }
    }

    // top front right (x=nNodes0-1,y=0,z=nNodes2-1)
    x = nNodes0-1;
    y = 0;
    z = nNodes2-1;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
        }
      }
    }

    // top back left (x=0,y=nNodes1-1,z=nNodes2-1)
    x = 0;
    y = nNodes1-1;
    z = nNodes2-1;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
        }
      }
    }

    // top back right (x=nNodes0-1,y=nNodes1-1,z=nNodes2-1)
    x = nNodes0-1;
    y = nNodes1-1;
    z = nNodes2-1;
    dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilCorner[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          massMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
        }
      }
    }

    // call MatAssemblyBegin, MatAssemblyEnd
    massMatrix->assembly(MAT_FINAL_ASSEMBLY);
  }
}

};