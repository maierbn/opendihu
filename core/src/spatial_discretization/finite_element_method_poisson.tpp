#include "spatial_discretization/finite_element_method.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "control/python_utility.h"


namespace SpatialDiscretization
{
  
template<class MeshType, class BasisFunctionType>
FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
FiniteElementMethodBaseRhs(DihuContext &context)
  : FiniteElementMethodBase<MeshType, BasisFunctionType>(context), timeSteppingInitialized_(false)
{
}

template<typename MeshType, typename BasisFunctionType, typename Term>
FiniteElementMethod<MeshType, BasisFunctionType, Term, Equation::hasLaplaceOperatorWithRhs<Term>>::
FiniteElementMethod(DihuContext &context) 
  : FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>(context)
{
}


template<class Mesh, class BasisFunction>
void FiniteElementMethodBaseRhs<Mesh, BasisFunction>::
setRightHandSide()
{
  LOG(DEBUG)<<"setRightHandSide for Poisson equation";

  int nDegreesOfFreedom = this->data_.mesh()->nDegreesOfFreedom();
  
  // fill rhs vector
  PetscErrorCode ierr;
 
  Vec &rightHandSide = this->data_.rightHandSide();

  // loop over rhs list
 
  // get the first dirichlet boundary condition from the list
  double rhsValue = PythonUtility::getOptionListBegin<double>(this->specificSettings_, "rightHandSide");
  node_idx_t nodeNo = 0; 
 
  // loop over Dirichlet boundary conditions and store them in rhs
  for (;
       !PythonUtility::getOptionListEnd(this->specificSettings_, "rightHandSide")
       && nodeNo < nDegreesOfFreedom; 
       PythonUtility::getOptionListNext<double>(this->specificSettings_, "rightHandSide", rhsValue), nodeNo++)
  {
    //                 vector         row     value
    ierr = VecSetValue(rightHandSide, nodeNo, rhsValue, INSERT_VALUES); CHKERRV(ierr);
  
  }
 
  // set remaining values to 0
  for (; nodeNo < nDegreesOfFreedom; nodeNo++)
  {
    //                 vector         row     value
    ierr = VecSetValue(rightHandSide, nodeNo, 0.0, INSERT_VALUES); CHKERRV(ierr);
  }
  this->multiplyRhsFactor();
}


template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
multiplyRhsFactor()
{
  LOG(FATAL) << "not implemented";
}

// 1D rhs
template<>
void FiniteElementMethodBaseRhs<Mesh::RegularFixed<1>, BasisFunction::Lagrange>::
multiplyRhsFactor()
{
  LOG(DEBUG)<<"multiplyRhsFactor (1D)";
 
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
  const double stencilCenter[3] = {1./6.*1, 1./6.*4, 1./6.*1};
  const double stencilSide[2] = {1./6.*2, 1./6.*1};
  
  // get all values 
  int nEntries;
  VecGetSize(rightHandSide, &nEntries);
  
  std::vector<int> indices(nEntries);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<double> vectorValues(nEntries);
  
  VecGetValues(rightHandSide, nEntries, indices.data(), vectorValues.data());
   
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
}

// 2D rhs
template<>
void FiniteElementMethodBaseRhs<Mesh::RegularFixed<2>, BasisFunction::Lagrange>::
multiplyRhsFactor()
{
  LOG(DEBUG)<<"multiplyRhsFactor (2D)";
 
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
  
  std::vector<int> indices(nEntries);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<double> vectorValues(nEntries);
  
  VecGetValues(rightHandSide, nEntries, indices.data(), vectorValues.data());
   
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
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
createDiscretizationMatrix()
{
}

template<>
void FiniteElementMethodBaseRhs<Mesh::RegularFixed<1>, BasisFunction::Lagrange>::
createDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(DEBUG)<<"createDiscretizationMatrix 1D";
 
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
void FiniteElementMethodBaseRhs<Mesh::RegularFixed<2>, BasisFunction::Lagrange>::
createDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(DEBUG)<<"createDiscretizationMatrix 2D";
      
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
 

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
initialize()
{
  FiniteElementMethodBase<MeshType, BasisFunctionType>::initialize(); 
}

template<class MeshType, class BasisFunctionType>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
mesh()
{
  return FiniteElementMethodBase<MeshType, BasisFunctionType>::mesh();
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
recoverRightHandSide(Vec &result)
{
  // create the discretization matrix if it does not already exist
  this->createDiscretizationMatrix();
  
  // dmatrix * f_strong = rhs_weak
  Vec &rhs = this->data_.rightHandSide();   // rhs in weak formulation
  Mat &dmatrix = this->data_.discretizationMatrix();
  
  LOG(DEBUG) << "recoverRightHandSide";
  
  PetscErrorCode ierr;
  
  // create linear solver context
  KSP ksp; 
  ierr = KSPCreate (PETSC_COMM_WORLD, &ksp); CHKERRV(ierr);  
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (ksp, dmatrix, dmatrix); CHKERRV(ierr);
  
  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC (ksp, &pc); CHKERRV(ierr);
  
  // set preconditioner type
  ierr = PCSetType (pc, PCJACOBI); CHKERRV(ierr);
  
  // set solution tolerances
  double relativeTolerance = PythonUtility::getOptionDouble(this->specificSettings_, "relativeTolerance", 1e-5, PythonUtility::Positive);

  //                            relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, relativeTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);

  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(result, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRV(ierr);
  
  // solve the system
  ierr = KSPSolve(ksp, rhs, result); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(ksp, &residualNorm); CHKERRV(ierr);
  
  LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output)
{
  if (!this->timeSteppingInitialized_)
  {
    this->setStiffnessMatrix();
    this->data_.finalAssembly();
    
    this->applyBoundaryConditions();
    
    this->timeSteppingInitialized_ = true;
  }
  
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  Vec &rhs = this->data_.rightHandSide();
  
  // compute rhs = stiffnessMatrix*input
  MatMult(stiffnessMatrix, input, rhs);
  
  recoverRightHandSide(output);
  
  this->data_.print();  
}

};