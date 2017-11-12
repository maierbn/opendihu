#include "spatial_discretization/finite_element_method.h"

#include <Python.h>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "mesh/regular_fixed.h"
#include "basis_function/lagrange.h"

#include <Python.h>
#include <memory>

#include "spatial_discretization/spatial_discretization.h"
#include "time_stepping_scheme/discretizable_in_time.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "control/petsc_utility.h"
#include "data_management/finite_elements.h"
#include "equation/laplace.h"
#include "equation/poisson.h"
#include "equation/type_traits.h"
#include "mesh/mesh.h"



namespace SpatialDiscretization
{
  
// 1D stiffness matrix
template<>
void FiniteElementMethodBase<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange>::
setStiffnessMatrix()
{
  LOG(DEBUG)<<"setStiffnessMatrix 1D";
 
  // get settings values
  int nElements = std::static_pointer_cast<Mesh::RegularFixed<1>>(data_.mesh())->nElements();
  double elementLength = std::static_pointer_cast<Mesh::RegularFixed<1>>(data_.mesh())->meshWidth(0);
  
  int nDegreesOfFreedom = data_.mesh()->nDegreesOfFreedom();
  
  LOG(DEBUG) << "Use settings nElements="<<nElements<<", elementLength="<<elementLength;
 
  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * (dxi/ds)^2 ds = l
  PetscErrorCode ierr;
 
  Mat &stiffnessMatrix = data_.stiffnessMatrix();
    
  // stencil values
  // stencil in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])
  const int center = 1;
  const double stencilCenter[3] = {1.0, -2.0, 1.0};
  const double stencilSide[2] = {-1.0, 1.0};
  
  for (node_idx_t dofNo = 0; dofNo < nDegreesOfFreedom; dofNo++)
  {
    // stencil for -Δu in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])
   
    //                 matrix           row        column
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo, stencilCenter[center]*elementLength, INSERT_VALUES); CHKERRV(ierr);
   
    if (dofNo+1 < nDegreesOfFreedom)
    {
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo+1, stencilCenter[center+1]*elementLength, INSERT_VALUES); CHKERRV(ierr);
    }
    if (dofNo-1 >= 0)
    {
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo-1, stencilCenter[center-1]*elementLength, INSERT_VALUES); CHKERRV(ierr);
    }
  }
  
  // set center values for boundaries
  ierr = MatSetValue(stiffnessMatrix, 0, 0, stencilSide[0]*elementLength, INSERT_VALUES); CHKERRV(ierr);
  ierr = MatSetValue(stiffnessMatrix, nDegreesOfFreedom-1, nDegreesOfFreedom-1, 
                     stencilSide[0]*elementLength, INSERT_VALUES); CHKERRV(ierr);
}

// 2D stiffness matrix
template<>
void FiniteElementMethodBase<Mesh::RegularFixed<2>, BasisFunction::Lagrange>::
setStiffnessMatrix()
{
  LOG(DEBUG)<<"setStiffnessMatrix 2D";
 
  // get settings values
  int nElements0 = std::static_pointer_cast<Mesh::RegularFixed<1>>(data_.mesh())->nElements(0);
  int nElements1 = std::static_pointer_cast<Mesh::RegularFixed<1>>(data_.mesh())->nElements(1);
  int nNodes0 = nElements0 + 1;
  int nNodes1 = nElements1 + 1;
  double elementLength0 = std::static_pointer_cast<Mesh::RegularFixed<1>>(data_.mesh())->meshWidth(0);
  double elementLength1 = std::static_pointer_cast<Mesh::RegularFixed<1>>(data_.mesh())->meshWidth(1);
  double integralFactor = elementLength0*elementLength1;
  
  //int nDegreesOfFreedom = data_.mesh()->nDegreesOfFreedom();
  
  LOG(DEBUG) << "Use settings nElements="<<nElements0<<"x"<<nElements1<<", elementLength="<<elementLength0<<"x"<<elementLength1;
  LOG(DEBUG) << "integralFactor="<<integralFactor;
  
  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * (dxi/ds)^2 ds = l

  // stencil for -Δu in 2D:     [1  1   1] (element contribution: [  1/6  1/3])
  //                        1/3*[1 _-8_ 1]                        [_-2/3_ 1/6]
  //                            [1  1   1]
  
  PetscErrorCode ierr;
  Mat &stiffnessMatrix = data_.stiffnessMatrix();
  
  const int center = 1;
  const double stencilCenter[3][3] = {
    {1./3, 1./3, 1./3},
    {1./3, -8./3, 1./3},
    {1./3, 1./3, 1./3}};
 
  auto dofIndex = [&nElements0, &nElements1](int x, int y){return y*(nElements0+1) + x;};
  double value = 0.0;
  
  // set entries for interior nodes
  for (int y=1; y<nNodes1-1; y++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      node_idx_t dofNo = dofIndex(x,y);
    
      // diagonal entry
      value = stencilCenter[center][center]*integralFactor;
      //                 matrix           row    column
      ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
     
      // left
      if (x > 0)
      {
        value = stencilCenter[center][center-1]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo+1, value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // right
      if (x < nNodes0-1)
      {
        value = stencilCenter[center][center+1]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo-1, value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // bottom
      if (y > 0)
      {
        value = stencilCenter[center-1][center]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y-1), value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // top
      if (y < nNodes1-1)
      {
        value = stencilCenter[center+1][center]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y+1), value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // bottom left
      if (y > 0 && x > 0)
      {
        value = stencilCenter[center-1][center-1]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y-1), value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // bottom right
      if (y > 0 && x < nNodes0-1)
      {
        value = stencilCenter[center-1][center+1]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y-1), value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // top left
      if (y < nNodes1-1 && x > 0)
      {
        value = stencilCenter[center+1][center-1]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y+1), value, INSERT_VALUES); CHKERRV(ierr);
      }
      
      // top right
      if (y < nNodes1-1 && x < nNodes0-1)
      {
        value = stencilCenter[center+1][center+1]*integralFactor;
        ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y+1), value, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }
  
  // set entries for boundary nodes on edges
  // left boundary
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    node_idx_t dofNo = dofIndex(x,y);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo,             -4.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y+1),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y-1),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom 
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y),   1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // right
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y+1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top right
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y-1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom right
  }
  
  // right boundary
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    node_idx_t dofNo = dofIndex(x,y);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo,             -4.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y+1),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y-1),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom 
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y),   1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // left
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y+1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top left
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y-1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom left
  }
  
  // bottom boundary
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    node_idx_t dofNo = dofIndex(x,y);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo,             -4.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // left
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // right
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y+1),   1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y+1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top left
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y+1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // top right
  }
  
  // top boundary
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    node_idx_t dofNo = dofIndex(x,y);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofNo,             -4.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // left
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y),   1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // right
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y-1),   1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y-1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom left
    ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y-1), 1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr); // bottom right
  } 
  
  // corner nodes
  // bottom left
  int x = 0;
  int y = 0;
  node_idx_t dofNo = dofIndex(x,y);
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y),     -2.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // self
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // right
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y+1),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // top
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y+1),  1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // top right
  
  // bottom right
  x = nNodes0-1;
  y = 0;
  dofNo = dofIndex(x,y);
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y),     -2.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // self
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // left
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y+1),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // top
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y+1),  1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // top left
  
  // top left
  x = 0;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y),     -2.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // self
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // right
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y-1),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // bottom
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x+1,y-1),  1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // bottom right
  
  // top right
  x = nNodes0-1;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y),     -2.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // self
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // left
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x,y-1),    1.0/6.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // bottom
  ierr = MatSetValue(stiffnessMatrix, dofNo, dofIndex(x-1,y-1),  1.0/3.*integralFactor, INSERT_VALUES); CHKERRV(ierr);    // bottom left
}
  
// 1D rhs
template<>
void FiniteElementMethodBaseRhs<Mesh::RegularFixed<1>, BasisFunction::Lagrange>::
transferRhsToWeakForm()
{
  LOG(DEBUG)<<"transferRhsToWeakForm (1D)";
 
  // get settings values
  int nElements = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->nElements();
  double elementLength = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->meshWidth(0);
  
  int nDegreesOfFreedom = this->data_.mesh()->nDegreesOfFreedom();
  
  LOG(DEBUG) << "Use settings nElements="<<nElements<<", elementLength="<<elementLength;
 
  // multiply factor to rhs
  // rhs *= stencil * elementLength
  PetscErrorCode ierr;
 
  Vec &rightHandSide = this->data_.rightHandSide();
    
  // stencil values
  // stencil in 1D: 1/6*[1 _4_ 1] (element contribution: 1/6*[_2_ 1])
  const int center = 1;
  const double stencilCenter[3] = {1./6., 4./6., 1./6.};
  const double stencilSide[2] = {2./6., 1./6.};
  
  // get all values 
  int nEntries;
  VecGetSize(rightHandSide, &nEntries);
  
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(rightHandSide, vectorValues);
   
  // loop over all dofs and set values with stencilCenter
  for (node_idx_t dofNo = 1; dofNo < nDegreesOfFreedom-1; dofNo++)
  {
    double value = 
      (stencilCenter[center-1]*vectorValues[dofNo-1] 
      + stencilCenter[center]*vectorValues[dofNo] 
      + stencilCenter[center+1]*vectorValues[dofNo+1]) * elementLength;
      
    ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  }
  
  // set values for boundaries with stencilSide
  node_idx_t dofNo = 0;
  double value = 
    (stencilSide[0]*vectorValues[dofNo] 
    + stencilSide[1]*vectorValues[dofNo+1]) * elementLength;
  ierr = VecSetValue(rightHandSide, 0, value, INSERT_VALUES); CHKERRV(ierr);
  
  dofNo = nDegreesOfFreedom-1;
  value = 
    (stencilSide[0]*vectorValues[dofNo]
    + stencilSide[1]*vectorValues[dofNo-1]) * elementLength;
  ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  
  VecAssemblyBegin(rightHandSide);
  VecAssemblyEnd(rightHandSide);
}

// 2D rhs
template<>
void FiniteElementMethodBaseRhs<Mesh::RegularFixed<2>, BasisFunction::Lagrange>::
transferRhsToWeakForm()
{
  LOG(DEBUG)<<"transferRhsToWeakForm (2D)";
 
  // get settings values
  int nElements0 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->nElements(0);
  int nElements1 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->nElements(1);
  int nNodes0 = nElements0 + 1;
  int nNodes1 = nElements1 + 1;
  double elementLength0 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->meshWidth(0);
  double elementLength1 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->meshWidth(1);
  double integralFactor = elementLength0*elementLength1;
  
  //int nDegreesOfFreedom = data_.mesh()->nDegreesOfFreedom();
  
  // multiply factor to rhs
  // rhs *= stencil * elementLength
  PetscErrorCode ierr;
 
  Vec &rightHandSide = this->data_.rightHandSide();
    
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
  double value;
 
  // get all values 
  int nEntries;
  VecGetSize(rightHandSide, &nEntries);
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(rightHandSide, vectorValues);
  
  // loop over all dofs and set values with stencilCenter
  // set entries for interior nodes
  for (int y=1; y<nNodes1-1; y++)
  {
    for (int x=1; x<nNodes0-1; x++)
    {
      value = 0;
      for (int i=-1; i<=1; i++)
      {
        for (int j=-1; j<=1; j++)
        {
          value += stencilCenter[center+i][center+j] * vectorValues[dofIndex(x+j, y+i)];
        }
      }
      value *= integralFactor;
        
      ierr = VecSetValue(rightHandSide, dofIndex(x,y), value, INSERT_VALUES); CHKERRV(ierr);
    }
  }
  
  // set entries for boundary nodes on edges
  // left boundary
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = 0;
    node_idx_t dofNo = dofIndex(x,y);
    
    value = 0;
    for (int i=-1; i<=1; i++)
    {
      value += stencilEdge[0][center+i] * vectorValues[dofIndex(x,y+i)]
        + stencilEdge[1][center+i] * vectorValues[dofIndex(x+1,y+i)];
    }
    value *= integralFactor;
    
    ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  }
  
  // right boundary
  for (int y=1; y<nNodes1-1; y++)
  {
    int x = nNodes0-1;
    node_idx_t dofNo = dofIndex(x,y);
    
    value = 0;
    for (int i=-1; i<=1; i++)
    {
      value += stencilEdge[0][center+i] * vectorValues[dofIndex(x,y+i)]
        + stencilEdge[1][center+i] * vectorValues[dofIndex(x-1,y+i)];
    }
    value *= integralFactor;
    
    ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  }
  
  // bottom boundary
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = 0;
    node_idx_t dofNo = dofIndex(x,y);
    
    value = 0;
    for (int i=-1; i<=1; i++)
    {
      value += stencilEdge[0][center+i] * vectorValues[dofIndex(x+i,y)]
        + stencilEdge[1][center+i] * vectorValues[dofIndex(x+i,y+1)];
    }
    value *= integralFactor;
    
    ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  }
  
  // top boundary
  for (int x=1; x<nNodes0-1; x++)
  {
    int y = nNodes1-1;
    node_idx_t dofNo = dofIndex(x,y);
    
    value = 0;
    for (int i=-1; i<=1; i++)
    {
      value += stencilEdge[0][center+i] * vectorValues[dofIndex(x+i,y)]
        + stencilEdge[1][center+i] * vectorValues[dofIndex(x+i,y-1)];
    }
    value *= integralFactor;
    
    ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  } 
 
  // corner nodes
  // bottom left
  int x = 0;
  int y = 0;
  node_idx_t dofNo = dofIndex(x,y);
  
  value = stencilCorner[0][0] * vectorValues[dofIndex(x,y)]
    + stencilCorner[0][1] * vectorValues[dofIndex(x+1,y)]
    + stencilCorner[1][0] * vectorValues[dofIndex(x,y+1)]
    + stencilCorner[1][1] * vectorValues[dofIndex(x+1,y+1)];
  value *= integralFactor;
  
  ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  
  // bottom right
  x = nNodes0-1;
  y = 0;
  dofNo = dofIndex(x,y);
  
  value = stencilCorner[0][0] * vectorValues[dofIndex(x,y)]
    + stencilCorner[0][1] * vectorValues[dofIndex(x-1,y)]
    + stencilCorner[1][0] * vectorValues[dofIndex(x,y+1)]
    + stencilCorner[1][1] * vectorValues[dofIndex(x-1,y+1)];
  value *= integralFactor;
  
  ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  
  // top left
  x = 0;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);
  
  value = stencilCorner[0][0] * vectorValues[dofIndex(x,y)]
    + stencilCorner[0][1] * vectorValues[dofIndex(x+1,y)]
    + stencilCorner[1][0] * vectorValues[dofIndex(x,y-1)]
    + stencilCorner[1][1] * vectorValues[dofIndex(x+1,y-1)];
  value *= integralFactor;
  
  ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  
  // top right
  x = nNodes0-1;
  y = nNodes1-1;
  dofNo = dofIndex(x,y);
  
  value = stencilCorner[0][0] * vectorValues[dofIndex(x,y)]
    + stencilCorner[0][1] * vectorValues[dofIndex(x-1,y)]
    + stencilCorner[1][0] * vectorValues[dofIndex(x,y-1)]
    + stencilCorner[1][1] * vectorValues[dofIndex(x-1,y-1)];
  value *= integralFactor;
  
  ierr = VecSetValue(rightHandSide, dofNo, value, INSERT_VALUES); CHKERRV(ierr);
  
  VecAssemblyBegin(rightHandSide);
  VecAssemblyEnd(rightHandSide);
}


template<>
void FiniteElementMethodBaseTimeStepping<Mesh::RegularFixed<1>, BasisFunction::Lagrange>::
createRhsDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(DEBUG)<<"createRhsDiscretizationMatrix 1D";
 
    // get settings values
    int nElements = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->nElements();
    double elementLength = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->meshWidth(0);
    
    int nDegreesOfFreedom = this->data_.mesh()->nDegreesOfFreedom();
    
    LOG(DEBUG) << "Use settings nElements="<<nElements<<", elementLength="<<elementLength;
   
    // multiply factor to rhs
    // rhs *= stencil * elementLength
    PetscErrorCode ierr;
   
    Mat &dmatrix = this->data_.discretizationMatrix();
      
    // dmatrix * f_strong = rhs_weak
    // row of dmatrix: contributions to a single entry in rhs_weak
    
    // stencil values
    // stencil in 1D: 1/6*[1 _4_ 1] (element contribution: 1/6*[_2_ 1])
    const int center = 1;
    const double stencilCenter[3] = {1./6.*1, 1./6.*4, 1./6.*1};
    const double stencilSide[2] = {1./6.*2, 1./6.*1};
    
    // loop over all dofs and set values in dmatrix with stencilCenter
    for (node_idx_t dofNo = 1; dofNo < nDegreesOfFreedom-1; dofNo++)
    { 
      ierr = MatSetValue(dmatrix, dofNo, dofNo-1, stencilCenter[center-1] * elementLength, INSERT_VALUES); CHKERRV(ierr);
      ierr = MatSetValue(dmatrix, dofNo, dofNo,   stencilCenter[center]   * elementLength, INSERT_VALUES); CHKERRV(ierr);
      ierr = MatSetValue(dmatrix, dofNo, dofNo+1, stencilCenter[center+1] * elementLength, INSERT_VALUES); CHKERRV(ierr);
    }
    
    // set values for boundaries with stencilSide
    node_idx_t dofNo = 0;    
    ierr = MatSetValue(dmatrix, dofNo, dofNo,   stencilSide[0] * elementLength, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofNo+1, stencilSide[1] * elementLength, INSERT_VALUES); CHKERRV(ierr);
    
    dofNo = nDegreesOfFreedom-1;
    ierr = MatSetValue(dmatrix, dofNo, dofNo,   stencilSide[0] * elementLength, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofNo-1, stencilSide[1] * elementLength, INSERT_VALUES); CHKERRV(ierr);
    
    ierr = MatAssemblyBegin(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  }
}

template<>
void FiniteElementMethodBaseTimeStepping<Mesh::RegularFixed<2>, BasisFunction::Lagrange>::
createRhsDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(DEBUG)<<"createRhsDiscretizationMatrix 2D";
      
    // get settings values
    int nElements0 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->nElements(0);
    int nElements1 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->nElements(1);
    int nNodes0 = nElements0 + 1;
    int nNodes1 = nElements1 + 1;
    double elementLength0 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->meshWidth(0);
    double elementLength1 = std::static_pointer_cast<Mesh::RegularFixed<1>>(this->data_.mesh())->meshWidth(1);
    double integralFactor = elementLength0*elementLength1;
      
    // multiply factor to rhs
    // rhs *= stencil * elementLength
    PetscErrorCode ierr;
   
    Mat &dmatrix = this->data_.discretizationMatrix();
      
    // dmatrix * f_strong = rhs_weak
    // row of dmatrix: contributions to a single entry in rhs_weak
    
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
        node_idx_t dofNo = dofIndex(x,y);
        for (int i=-1; i<=1; i++)
        {
          for (int j=-1; j<=1; j++)
          {
            node_idx_t secondDofNo = dofIndex(x+j, y+i);
            double factor = stencilCenter[center+i][center+j] * integralFactor;
            ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
          }
        }
      }
    }
    
    // set entries for boundary nodes on edges
    // left boundary
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = 0;
      node_idx_t dofNo = dofIndex(x,y);
      
      for (int i=-1; i<=1; i++)
      { 
        node_idx_t secondDofNo = dofIndex(x, y+i);
        double factor = stencilEdge[0][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
        
        secondDofNo = dofIndex(x+1, y+i);
        factor = stencilEdge[1][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
      }
    }
    
    // right boundary
    for (int y=1; y<nNodes1-1; y++)
    {
      int x = nNodes0-1;
      node_idx_t dofNo = dofIndex(x,y);
      
      for (int i=-1; i<=1; i++)
      {
        node_idx_t secondDofNo = dofIndex(x, y+i);
        double factor = stencilEdge[0][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
        
        secondDofNo = dofIndex(x-1, y+i);
        factor = stencilEdge[1][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
      }
    }
    
    // bottom boundary
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = 0;
      node_idx_t dofNo = dofIndex(x,y);
      
      for (int i=-1; i<=1; i++)
      {
        node_idx_t secondDofNo = dofIndex(x+i, y);
        double factor = stencilEdge[0][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
        
        secondDofNo = dofIndex(x+i, y+1);
        factor = stencilEdge[1][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
      }
    }
    
    // top boundary
    for (int x=1; x<nNodes0-1; x++)
    {
      int y = nNodes1-1;
      node_idx_t dofNo = dofIndex(x,y);
      
      for (int i=-1; i<=1; i++)
      {
        node_idx_t secondDofNo = dofIndex(x+i, y);
        double factor = stencilEdge[0][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
        
        secondDofNo = dofIndex(x+i, y-1);
        factor = stencilEdge[1][center+i] * integralFactor;
        ierr = MatSetValue(dmatrix, dofNo, secondDofNo, factor, INSERT_VALUES); CHKERRV(ierr);
      }
    } 
   
    // corner nodes
    // bottom left
    int x = 0;
    int y = 0;
    node_idx_t dofNo = dofIndex(x,y);
    
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x+1,y),   stencilCorner[0][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y+1),   stencilCorner[1][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x+1,y+1), stencilCorner[1][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
        
    // bottom right
    x = nNodes0-1;
    y = 0;
    dofNo = dofIndex(x,y);
    
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x-1,y),   stencilCorner[0][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y+1),   stencilCorner[1][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x-1,y+1), stencilCorner[1][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
        
    // top left
    x = 0;
    y = nNodes1-1;
    dofNo = dofIndex(x,y);
    
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x+1,y),   stencilCorner[0][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y-1),   stencilCorner[1][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x+1,y-1), stencilCorner[1][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
        
    // top right
    x = nNodes0-1;
    y = nNodes1-1;
    dofNo = dofIndex(x,y);
    
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y),     stencilCorner[0][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x-1,y),   stencilCorner[0][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x,y-1),   stencilCorner[1][0] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
    ierr = MatSetValue(dmatrix, dofNo, dofIndex(x-1,y-1), stencilCorner[1][1] * integralFactor, INSERT_VALUES); CHKERRV(ierr);
  
   
    ierr = MatAssemblyBegin(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    
  }
}
 


};    // namespace