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
//#include "time_stepping_scheme/discretizable_in_time.h"
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
  element_no_t nElements = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->nLocalElements();
  double elementLength = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->meshWidth();
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  double factor = prefactor*1./elementLength;

  dof_no_t nDegreesOfFreedom = this->data_.mesh()->nLocalNodes();

  LOG(DEBUG) << "  Use settings nElements="<<nElements<<", elementLength="<<elementLength;

  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * (dxi/ds)^2 ds = l
  PetscErrorCode ierr;

  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  // stencil values
  // stencil in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])
  const int center = 1;
  const double stencilCenter[3] = {1.0, -2.0, 1.0};
  const double stencilSide[2] = {-1.0, 1.0};

  for (node_no_t dofNo = 0; dofNo < nDegreesOfFreedom; dofNo++)
  {
    // stencil for -Δu in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])

    //                 matrix           row        column
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo, stencilCenter[center]*factor, INSERT_VALUES); CHKERRV(ierr);

    if (dofNo+1 < nDegreesOfFreedom)
    {
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo+1, stencilCenter[center+1]*factor, INSERT_VALUES); CHKERRV(ierr);
    }
    if (dofNo-1 >= 0)
    {
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo-1, stencilCenter[center-1]*factor, INSERT_VALUES); CHKERRV(ierr);
    }
  }

  // set center values for boundaries
  ierr = MatSetValue(stiffnessMatrix, 0, 0, stencilSide[0]*factor, INSERT_VALUES); CHKERRV(ierr);
  ierr = MatSetValue(stiffnessMatrix, nDegreesOfFreedom-1, nDegreesOfFreedom-1,
                     stencilSide[0]*factor, INSERT_VALUES); CHKERRV(ierr);
}

// 2D stiffness matrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<2>, Equation::hasLaplaceOperator<Term>,BasisFunction::LagrangeOfOrder<1>>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix 2D for Mesh::RegularFixed";

  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  // get settings value
  element_no_t nElements0 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->nElementsPerCoordinateDirectionLocal(0);
  element_no_t nElements1 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->nElementsPerCoordinateDirectionLocal(1);
  node_no_t nNodes0 = nElements0 + 1;
  node_no_t nNodes1 = nElements1 + 1;
  double elementLength0 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->meshWidth();
  double elementLength1 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->meshWidth();
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

  PetscErrorCode ierr;
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
  node_no_t dofNo;

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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j), value, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }
  }

  // set entries for boundary nodes on edges
  // left boundary (x=0)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    node_no_t dofNo = dofIndex(x,y);

    for (int i=-1; i<=0; i++) // -x
    {
      for (int j=-1; j<=1; j++) // y
      {
        value = stencilEdge[center+i][center+j]*integralFactor;
        //                 matrix           row    column
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j), value, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }

  // right boundary (x=nNodes0-1)
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    node_no_t dofNo = dofIndex(x,y);
    for (int i=-1; i<=0; i++) // x
    {
      for (int j=-1; j<=1; j++) // y
      {
        value = stencilEdge[center+i][center+j]*integralFactor;
        //                 matrix           row    column
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j), value, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }

  // bottom boundary (y=0)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    node_no_t dofNo = dofIndex(x,y);
    for (int i=-1; i<=1; i++) // x
    {
      for (int j=-1; j<=0; j++) // -y
      {
        value = stencilEdge[center+j][center+i]*integralFactor;
        //                 matrix           row    column
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j), value, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }

  // top boundary (y=nNodes1-1)
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    node_no_t dofNo = dofIndex(x,y);
    for (int i=-1; i<=1; i++) // x
    {
      for (int j=-1; j<=0; j++) // y
      {
        value = stencilEdge[center+j][center+i]*integralFactor;
        //                 matrix           row    column
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j), value, INSERT_VALUES); CHKERRV(ierr);
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
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y-j), value, INSERT_VALUES); CHKERRV(ierr);
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
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j), value, INSERT_VALUES); CHKERRV(ierr);
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
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j), value, INSERT_VALUES); CHKERRV(ierr);
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
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j), value, INSERT_VALUES); CHKERRV(ierr);
    }
  }
}

// 3D stiffness matrix
template<typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<3>, Equation::hasLaplaceOperator<Term>,BasisFunction::LagrangeOfOrder<1>>::
setStiffnessMatrix()
{
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> BasisOnMeshType;

  LOG(TRACE)<<"setStiffnessMatrix 3D for Mesh::RegularFixed";

  // get settings values
  element_no_t nElements0 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->nElementsPerCoordinateDirectionLocal(0);
  element_no_t nElements1 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->nElementsPerCoordinateDirectionLocal(1);
  element_no_t nElements2 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->nElementsPerCoordinateDirectionLocal(2);
  node_no_t nNodes0 = nElements0 + 1;
  node_no_t nNodes1 = nElements1 + 1;
  node_no_t nNodes2 = nElements2 + 1;
  double elementLength0 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->meshWidth();
  double elementLength1 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->meshWidth();
  double elementLength2 = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh())->meshWidth();

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

  PetscErrorCode ierr;
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

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
              ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
            }
          }
        }
      }
    }
  }

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
            ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
            ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
            ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
            ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
            ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
            ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y-j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
          ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y-j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z-k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y-j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y-j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
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
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+i, y+j, z+k), value, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }
}



};    // namespace