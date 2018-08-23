#include "spatial_discretization/finite_element_method/02_stiffness_matrix.h"

#include <Python.h>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "mesh/structured_regular_fixed.h"
#include "basis_function/lagrange.h"

#include <Python.h>
#include <memory>

#include "spatial_discretization/spatial_discretization.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "data_management/finite_elements.h"
#include "equation/laplace.h"
#include "equation/poisson.h"
#include "equation/type_traits.h"
#include "mesh/mesh.h"
#include "control/types.h"



namespace SpatialDiscretization
{

// 1D stiffness matrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<1>, Equation::hasLaplaceOperator<Term>,BasisFunction::LagrangeOfOrder<1>>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix 1D for Mesh::RegularFixed";

  typedef typename BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // get settings values
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  element_no_t nElements = mesh->nElementsLocal();
  node_no_t nNodes0 = mesh->nNodesLocalWithGhosts(0);
  double elementLength = mesh->meshWidth();

  double integralFactor = 1./elementLength;
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  integralFactor = prefactor*integralFactor;

  LOG(DEBUG) << "  Use settings nElements="<<nElements<<", elementLength="<<elementLength;

  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * (dxi/ds)^2 ds = l
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  // stencil values
  // stencil for -Δu in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])
  const int center = 1;
  const double stencilCenter[3] = {1.0, -2.0, 1.0};
  const double stencilSide[2] = {1.0, -1.0};  // left side of stencil
  
  double value;
  dof_no_t dofNo;

  // loop over all dofs and set values with stencilCenter
  // set entries for interior nodes
  for (int x=1; x<nNodes0-1; x++)
  {
    dofNo = x;
    for (int i=-1; i<=1; i++) // x
    {
      value = stencilCenter[center+i]*integralFactor;
      //                 matrix           row    column
      stiffnessMatrix->setValue(dofNo, x+i, value, INSERT_VALUES);
    }
  }

  // call MatAssemblyBegin, MatAssemblyEnd
  stiffnessMatrix->assembly(MAT_FLUSH_ASSEMBLY);

  // set entries for boundary nodes on edges
  // left boundary (x=0)
  int x = 0;
  dofNo = x;

  for (int i=-1; i<=0; i++) // -x
  {
    value = stencilSide[center+i]*integralFactor;
    stiffnessMatrix->setValue(dofNo, x-i, value, ADD_VALUES);
  }

  // right boundary (x=nNodes0-1)
  x = nNodes0-1;
  dofNo = x;
  for (int i=-1; i<=0; i++) // x
  {
    value = stencilSide[center+i]*integralFactor;
    stiffnessMatrix->setValue(dofNo, x+i, value, ADD_VALUES);
  }
  
  // call MatAssemblyBegin, MatAssemblyEnd
  //stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);
}

// 2D stiffness matrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<2>, Equation::hasLaplaceOperator<Term>,BasisFunction::LagrangeOfOrder<1>>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix 2D for Mesh::RegularFixed";

  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // get settings value
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  element_no_t nElements0 = mesh->nElementsPerCoordinateDirectionLocal(0);
  element_no_t nElements1 = mesh->nElementsPerCoordinateDirectionLocal(1);
  node_no_t nNodes0 = mesh->nNodesLocalWithGhosts(0);
  node_no_t nNodes1 = mesh->nNodesLocalWithGhosts(1);
  double elementLength0 = mesh->meshWidth();
  double elementLength1 = mesh->meshWidth();
  if (fabs(elementLength0-elementLength1) > 1e-15)
  {
    LOG(ERROR) << "Mesh resolution of 2D regular fixed mesh is not uniform! " << std::endl
      << "Mesh widths: x: " << elementLength0 << ", y: " << elementLength1 << std::endl
      << "This means that the stiffness matrix will be wrong. Mesh::RegularFixed meshes use stencil notation "
      << "and can only handle uniform meshes correctly. To use non-uniform meshes, consider using Mesh::Deformable!";
  }

  double integralFactor = 1.;
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  integralFactor = prefactor*integralFactor;

  LOG(DEBUG) << "Use settings nElements="<<nElements0<<"x"<<nElements1<<", elementLength="<<elementLength0<<"x"<<elementLength1;
  LOG(DEBUG) << "integralFactor="<<integralFactor;

  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * (dxi/ds)^2 ds = l

  // stencil for -Δu in 2D:     [1  1   1] (element contribution: [  1/6  1/3])
  //                        1/3*[1 _-8_ 1]                        [_-2/3_ 1/6]
  //                            [1  1   1]

  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  const int center = 1;
  const double stencilCenter[3][3] = {
    {1./3, 1./3, 1./3},
    {1./3, -8./3, 1./3},
    {1./3, 1./3, 1./3}};

  const double stencilEdge[2][3] = {
    {1./3, 1./3, 1./3},
    {1./6, -4./3, 1./6}
  };

  const double stencilCorner[2][2] = {
    {1./3, 1./6},
    {1./6, -2./3}
  };

  auto dofIndex = [&nNodes0](int x, int y){return y*nNodes0 + x;};
  double value;
  dof_no_t dofNo;

  // loop over all dofs and set values with stencilCenter
  // set entries for interior nodes
  for (int y=1; y<nNodes1-1; y++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      dofNo = dofIndex(x, y);
      for (int i=-1; i<=1; i++) // x
      {
        for (int j=-1; j<=1; j++) // y
        {
          value = stencilCenter[center+i][center+j]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j), value, INSERT_VALUES);
        }
      }
    }
  }

  // call MatAssemblyBegin, MatAssemblyEnd
  stiffnessMatrix->assembly(MAT_FLUSH_ASSEMBLY);
  
  // set entries for boundary nodes on edges
  // left boundary (x=0)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    dof_no_t dofNo = dofIndex(x,y);

    for (int i=-1; i<=0; i++) // -x
    {
      for (int j=-1; j<=1; j++) // y
      {
        value = stencilEdge[center+i][center+j]*integralFactor;
        //                 matrix           row    column
        stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j), value, ADD_VALUES);
      }
    }
  }

  // right boundary (x=nNodes0-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    dof_no_t dofNo = dofIndex(x,y);
    for (int i=-1; i<=0; i++) // x
    {
      for (int j=-1; j<=1; j++) // y
      {
        value = stencilEdge[center+i][center+j]*integralFactor;
        //                 matrix           row    column
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j), value, ADD_VALUES);
      }
    }
  }

  // bottom boundary (y=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    dof_no_t dofNo = dofIndex(x,y);
    for (int i=-1; i<=1; i++) // x
    {
      for (int j=-1; j<=0; j++) // -y
      {
        value = stencilEdge[center+j][center+i]*integralFactor;
        //                 matrix           row    column
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j), value, ADD_VALUES);
      }
    }
  }

  // top boundary (y=nNodes1-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    dof_no_t dofNo = dofIndex(x,y);
    for (int i=-1; i<=1; i++) // x
    {
      for (int j=-1; j<=0; j++) // y
      {
        value = stencilEdge[center+j][center+i]*integralFactor;
        //                 matrix           row    column
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j), value, ADD_VALUES);
      }
    }
  }

  // corner nodes
  int x,y;

  // bottom left (x=0, y=0)
  x = 0;
  y = 0;
  dofNo = dofIndex(x,y);

  for (int i=-1; i<=0; i++) // -x
  {
    for (int j=-1; j<=0; j++) // -y
    {
      value = stencilCorner[center+i][center+j]*integralFactor;
      //                 matrix           row    column
      stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y-j), value, ADD_VALUES);
    }
  }

  // bottom right (x=nNodes0-1, y=0)
  x = nNodes0-1;
  y = 0;
  dofNo = dofIndex(x,y);

  for (int i=-1; i<=0; i++) // x
  {
    for (int j=-1; j<=0; j++) // -y
    {
      value = stencilCorner[center+i][center+j]*integralFactor;
      //                 matrix           row    column
      stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j), value, ADD_VALUES);
    }
  }

  // top left (x=0, y=nNodes1-1)
  x = 0;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);

  for (int i=-1; i<=0; i++) // -x
  {
    for (int j=-1; j<=0; j++) // y
    {
      value = stencilCorner[center+i][center+j]*integralFactor;
      //                 matrix           row    column
      stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j), value, ADD_VALUES);
    }
  }

  // top right (x=nNodes0-1, y=nNodes1-1)
  x = nNodes0-1;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);

  for (int i=-1; i<=0; i++) // x
  {
    for (int j=-1; j<=0; j++) // y
    {
      value = stencilCorner[center+i][center+j]*integralFactor;
      //                 matrix           row    column
      stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j), value, ADD_VALUES);
    }
  }
  
  //stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);
}

// 3D stiffness matrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<3>, Equation::hasLaplaceOperator<Term>,BasisFunction::LagrangeOfOrder<1>>::
setStiffnessMatrix()
{
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  LOG(TRACE)<<"setStiffnessMatrix 3D for Mesh::RegularFixed";

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

  if (fabs(elementLength0-elementLength1) > 1e-15 || fabs(elementLength0-elementLength2) > 1e-15)
  {
    LOG(ERROR) << "Mesh resolution of 3D regular fixed mesh is not uniform! " << std::endl
      << "Mesh widths: x: " << elementLength0 << ", y: " << elementLength1 << ", z: " << elementLength2 << std::endl
      << "This means that the stiffness matrix will be wrong. Mesh::RegularFixed meshes use stencil notation "
      << "and can only handle uniform meshes correctly. To use non-uniform meshes, consider using Mesh::Deformable!";
  }

  double integralFactor = elementLength0;
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  integralFactor = prefactor*integralFactor;

  LOG(DEBUG) << "Use settings nElements="<<nElements0<<"x"<<nElements1<<"x"<<nElements2<<
    ", elementLength="<<elementLength0<<"x"<<elementLength1<<"x"<<elementLength2;
  LOG(DEBUG) << "integralFactor="<<integralFactor;

  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * (dxi/ds)^2 ds = l

  // stencil for -Δu in 3D: bottom: [1  2  1]   (element contribution:  center: [    0  1/12]
  //                           1/12*[2 _0_ 2]                                   [_-4/12_   0]
  //                                [1  2  1]
  //                                                                     bottom:[ 1/12  1/12]
  //                        center: [ 2    0  2]                                [    0  1/12] )
  //                           1/12*[ 0 _-32_ 0]
  //                                [ 2    0  2]
  //
  //                        top: like bottom

  // coordinate system
  // x axis: left -> right
  // y axis: front -> back
  // z axis: bottom -> top


  const int center = 1;
  const double stencilCenter[3][3][3] = {
    {{1./12, 2./12, 1./12},   //bottom
    {2./12, 0./12, 2./12},
    {1./12, 2./12, 1./12}},
    {{2./12, 0./12, 2./12},   //center
    {0./12, -32./12, 0./12},
    {2./12, 0./12, 2./12}},
    {{1./12, 2./12, 1./12},   //top
    {2./12, 0./12, 2./12},
    {1./12, 2./12, 1./12}},
  };

  const double stencilBoundarySurface[2][3][3] = {
    {{1./12, 2./12, 1./12},   //bottom
    {2./12, 0./12, 2./12},
    {1./12, 2./12, 1./12}},
    {{1./12, 0./12, 1./12},   //center
    {0./12, -16./12, 0./12},
    {1./12, 0./12, 1./12}},
  };
  const double stencilBoundaryEdge[2][2][3] = {
    {{1./12, 2./12, 1./12},   //bottom
    {1./12, 0./12, 1./12}},
    {{1./12, 0./12, 1./12},
    {0./12, -8./12, 0./12}}    //center
  };

  const double stencilCorner[2][2][2] = {
    {{1./12, 1./12},
    {1./12, 0./12}},    //bottom
    {{1./12, 0./12},
    {0./12, -4./12}},    //center
  };

  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  auto dofIndex = [&nNodes0, &nNodes1](int x, int y, int z){
    return z*nNodes0*nNodes1 + y*nNodes0 + x;
  };
  double value;
  dof_no_t dofNo;

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
              stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES);
            }
          }
        }
      }
    }
  }

  // call MatAssemblyBegin, MatAssemblyEnd
  stiffnessMatrix->assembly(MAT_FLUSH_ASSEMBLY);
  
  // set entries for boundary nodes on surface boundaries
  // left boundary (x = 0)
  for (int z=1; z<nNodes2-1; z++)
  {
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = 0;
      dof_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // -x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundarySurface[center+i][center+j][center+k]*integralFactor;
            //                 matrix           row    column
            stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
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
      dof_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=0; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundarySurface[center+i][center+j][center+k]*integralFactor;
            //                 matrix           row    column
            stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
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
      dof_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // -y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundarySurface[center+j][center+i][center+k]*integralFactor;
            //                 matrix           row    column
            stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
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
      dof_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=0; j++)   // y
        {
          for (int k=-1; k<=1; k++)   // z
          {
            value = stencilBoundarySurface[center+j][center+i][center+k]*integralFactor;
            //                 matrix           row    column
            stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
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
      dof_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // -z
          {
            value = stencilBoundarySurface[center+k][center+i][center+j]*integralFactor;
            //                 matrix           row    column
            stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
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
      dof_no_t dofNo = dofIndex(x,y,z);
      for (int i=-1; i<=1; i++)    // x
      {
        for (int j=-1; j<=1; j++)   // y
        {
          for (int k=-1; k<=0; k++)   // z
          {
            value = stencilBoundarySurface[center+k][center+i][center+j]*integralFactor;
            //                 matrix           row    column
            stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
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
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j, z-k), value, ADD_VALUES);
        }
      }
    }
  }

  // bottom right (x=nNodes0-1,z=0)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    int z = 0;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
        }
      }
    }
  }

  // top left (x=0,z=nNodes2-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    int z = nNodes2-1;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // top right (x=nNodes0-1,z=nNodes2-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    int z = nNodes2-1;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=1; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilBoundaryEdge[center+i][center+k][center+j]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // bottom front (y=0,z=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    int z = 0;
    dof_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j, z-k), value, ADD_VALUES);
        }
      }
    }
  }

  // bottom back (y=nNodes1-1,z=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    int z = 0;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // -z
        {
          value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
        }
      }
    }
  }

  // top front (y=0,z=nNodes2-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    int z = nNodes2-1;
    dof_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // top back (y=nNodes1-1,z=nNodes2-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    int z = nNodes2-1;
    dof_no_t dofNo = dofIndex(x,y,z);

    value = 0;
    for (int i=-1; i<=1; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=0; k++)   // z
        {
          value = stencilBoundaryEdge[center+j][center+k][center+i]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // left front (x=0,y=0)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = 0;
    int y = 0;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y-j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // left back (x=0,y=nNodes1-1)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = 0;
    int y = nNodes1-1;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // -x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // right front (x=nNodes0-1,y=0)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = nNodes0-1;
    int y = 0;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // -y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
        }
      }
    }
  }

  // right back (x=nNodes0-1,y=nNodes1-1)
  for (int z=1; z<nNodes2-1; z++)
  {
    int x = nNodes0-1;
    int y = nNodes1-1;
    dof_no_t dofNo = dofIndex(x,y,z);
    for (int i=-1; i<=0; i++)    // x
    {
      for (int j=-1; j<=0; j++)   // y
      {
        for (int k=-1; k<=1; k++)   // z
        {
          value = stencilBoundaryEdge[center+i][center+j][center+k]*integralFactor;
          //                 matrix           row    column
          stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y-j, z-k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j, z-k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j, z-k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z-k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y-j, z+k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y-j, z+k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x-i, y+j, z+k), value, ADD_VALUES);
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
        stiffnessMatrix->setValue(dofNo, dofIndex(x+i, y+j, z+k), value, ADD_VALUES);
      }
    }
  }
  
  //stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);
}



};    // namespace