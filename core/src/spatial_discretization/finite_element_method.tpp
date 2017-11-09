#include "spatial_discretization/finite_element_method.h"

#include <iostream>
#include <petscksp.h>

#include <Python.h>

#include "control/types.h"
#include "control/python_utility.h"
#include "easylogging++.h"


namespace SpatialDiscretization
{

template<typename Mesh, typename BasisFunction>
FiniteElementMethodBase<Mesh, BasisFunction>::FiniteElementMethodBase(DihuContext &context) :
  context_(context), data_(context)
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::extractDict(topLevelSettings, "FiniteElementMethod");
  
  data_.setMesh(std::make_shared<Mesh>(specificSettings_));
  if(data_.mesh())
  {
    LOG(DEBUG) << "mesh is set";
  }
  else
  {
    LOG(DEBUG) << "mesh is not set";
  }
  
  //LOG(DEBUG) << "exit in finite_element_method.tpp:31";
  //exit(0);
}


template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBase<MeshType, BasisFunctionType>::
applyBoundaryConditions()
{
  // PETSc Mat object for stiffness matrix needs to be assembled for this
 
  LOG(DEBUG)<<"applyBoundaryConditions for Laplace";
 
  int nDegreesOfFreedom = this->data_.mesh()->nDegreesOfFreedom();
  
  Vec &rightHandSide = data_.rightHandSide();
  Mat &stiffnessMatrix = data_.stiffnessMatrix();
  PetscErrorCode ierr;
  
  // add Dirichlet boundary conditions
  
  // get the first dirichlet boundary condition from the list
  std::pair<node_idx_t, double> boundaryCondition 
    = PythonUtility::getOptionDictBegin<node_idx_t, double>(specificSettings_, "DirichletBoundaryCondition");
  
  // loop over Dirichlet boundary conditions
  for (; !PythonUtility::getOptionDictEnd(specificSettings_, "DirichletBoundaryCondition"); 
       PythonUtility::getOptionDictNext<node_idx_t, double>(specificSettings_, "DirichletBoundaryCondition", boundaryCondition))
  {
    node_idx_t boundaryConditionNodeIndex = boundaryCondition.first;
    double boundaryConditionValue = boundaryCondition.second;
    
    if (boundaryConditionNodeIndex < 0)
      continue;
    
    if (boundaryConditionNodeIndex > nDegreesOfFreedom) 
    {
      LOG(WARNING) << "Boundary condition specified for degree of freedom no. "<<boundaryConditionNodeIndex
       <<", but scenario has only "<<nDegreesOfFreedom<<" degrees of freedom.";
       continue;
    }
    
    // set rhs entry to prescribed value
    ierr = VecSetValue(rightHandSide, boundaryConditionNodeIndex, boundaryConditionValue, INSERT_VALUES); CHKERRV(ierr);

    LOG(DEBUG) << "  BC node " << boundaryConditionNodeIndex << " value " << boundaryConditionValue;
    
    // get the column number boundaryConditionNodeIndex of the stiffness matrix. It is needed for updating the rhs.
    std::vector<int> rowIndices(nDegreesOfFreedom);
    std::iota (rowIndices.begin(), rowIndices.end(), 0);    // fill with increasing numbers: 0,1,2,...
    std::vector<int> columnIndices = {boundaryConditionNodeIndex};
    
    std::vector<double> coefficients(nDegreesOfFreedom);
    
    ierr = MatGetValues(stiffnessMatrix, nDegreesOfFreedom, rowIndices.data(), 1, columnIndices.data(), coefficients.data());
        
    // set values of row and column of the DOF to zero and diagonal entry to 1
    ierr = MatZeroRowsColumns(stiffnessMatrix, 1, &boundaryConditionNodeIndex, 1.0, NULL, NULL);  CHKERRV(ierr);

    // update rhs
    for (node_idx_t rowNo = 0; rowNo < nDegreesOfFreedom; rowNo++)
    {
      if (rowNo == boundaryConditionNodeIndex)
       continue;
     
      // update rhs value to be f_new = f_old - m_{ij}*u_{i} where i is the index of the prescribed node,
      // m_{ij} is entry of stiffness matrix and u_{i} is the prescribed value
      double rhsSummand = -coefficients[rowNo] * boundaryConditionValue;
      ierr = VecSetValue(rightHandSide, rowNo, rhsSummand, ADD_VALUES); CHKERRV(ierr);
      
      LOG_IF(false,DEBUG) << "  in row " << rowNo << " add " << rhsSummand << " to rhs, coefficient: " << coefficients[rowNo];
    }
  }
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBase<MeshType, BasisFunctionType>::
setStiffnessMatrix()
{
  LOG(FATAL) << "Not implemented!";
}

// 1D stiffness matrix
template<>
void FiniteElementMethodBase<Mesh::RegularFixed<1>, BasisFunction::Lagrange>::
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
  
template<typename MeshT, typename BasisFunction>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodBase<MeshT, BasisFunction>::mesh()
{
  return data_.mesh();
}
  
template<typename Mesh, typename BasisFunction>
void FiniteElementMethodBase<Mesh, BasisFunction>::initialize()
{
  setStiffnessMatrix();
  setRightHandSide();
  data_.finalAssembly();
  
  applyBoundaryConditions();
}
  
template<typename Mesh, typename BasisFunction>
void FiniteElementMethodBase<Mesh, BasisFunction>::run()
{
  initialize();
  solve();
  data_.print();
  
  context_.writeOutput(data_, specificSettings_);
}

template<typename Mesh, typename BasisFunction>
void FiniteElementMethodBase<Mesh, BasisFunction>::solve()
{
  LOG(DEBUG) << "FiniteElementMethod::solve";
  
  PetscErrorCode ierr;
  
  Mat &stiffnessMatrix = data_.stiffnessMatrix();
  
  // create linear solver context
  KSP ksp; 
  ierr = KSPCreate (PETSC_COMM_WORLD, &ksp); CHKERRV(ierr);  
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (ksp, stiffnessMatrix, stiffnessMatrix); CHKERRV(ierr);
  
  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC (ksp, &pc); CHKERRV(ierr);
  
  // set preconditioner type
  ierr = PCSetType (pc, PCJACOBI); CHKERRV(ierr);
  
  // set solution tolerances
  double relativeTolerance = PythonUtility::getOptionDouble(specificSettings_, "relativeTolerance", 1e-5, PythonUtility::Positive);

  //                            relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, relativeTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);

  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(data_.solution(), scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRV(ierr);
  
  // solve the system
  ierr = KSPSolve(ksp, data_.rightHandSide(), data_.solution()); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(ksp, &residualNorm); CHKERRV(ierr);
  
  LOG(INFO) << "Solution done in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  
  // check if solution is correct
  if (false)
  {
    // get rhs and solution from PETSc
    int vectorSize = 0;
    VecGetSize(data_.solution(), &vectorSize);

    std::vector<int> indices(vectorSize);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> solution(vectorSize);
    std::vector<double> rhs(vectorSize);

    VecGetValues(data_.solution(), vectorSize, indices.data(), solution.data());
    VecGetValues(data_.rightHandSide(), vectorSize, indices.data(), rhs.data());
    
    // get stiffness matrix
    int nRows, nColumns;
    MatGetSize(data_.stiffnessMatrix(), &nRows, &nColumns);
    std::vector<int> rowIndices(nRows);
    std::iota(rowIndices.begin(), rowIndices.end(), 0);
    std::vector<int> columnIndices(nColumns);
    std::iota(columnIndices.begin(), columnIndices.end(), 0);
    std::vector<double> matrixValues(nRows*nColumns);
    
    std::vector<long int> nEntries = {nRows, nColumns};
    
    MatGetValues(data_.stiffnessMatrix(), nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
     
    std::vector<double> f(vectorSize);
    
    // compute f = matrix * solution
    
    for(int i=0; i<vectorSize; i++)
    {
      f[i] = 0.0;
      for(int j=0; j<vectorSize; j++)
      {
        f[i] += matrixValues[i*nColumns + j] * solution[j];
      }
    }
    
    // compute residual norm
    double res = 0.0;
    for(int i=0; i<vectorSize; i++)
    {
      res += (f[i] - rhs[i]) * (f[i] - rhs[i]);
      LOG(DEBUG) << i << ". solution="<<solution[i]<<", f="<<f[i]<<", rhs="<<rhs[i]<<", squared error: "<<(f[i] - rhs[i]) * (f[i] - rhs[i]);
    }
    
    LOG(DEBUG) << "res=" << res;
  }
  
  // clean up ksp solver context
  KSPDestroy(&ksp);
}
  
};