#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include <iostream>
#include <petscksp.h>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"

#include "mesh/regular_fixed.h"
#include "basis_function/lagrange.h"
#include "mesh/mesh_manager.h"

namespace SpatialDiscretization
{

template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::
FiniteElementMethodBase(const DihuContext &context) :
  context_(context), data_(context)
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "FiniteElementMethod");
  outputWriterManager_.initialize(specificSettings_);
  
  LOG(DEBUG) << "FiniteElementMethodBase::FiniteElementMethodBase querying meshManager for mesh, specificSettings_:";
  PythonUtility::printDict(specificSettings_);
  
  data_.setMesh(context_.meshManager()->mesh<MeshType>(specificSettings_));
  if(data_.mesh())
    LOG(DEBUG) << "FiniteElementMethodBase: mesh is set";
  else
    LOG(DEBUG) << "FiniteElementMethodBase: mesh is not set";
}


template<class MeshType, class BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::
applyBoundaryConditions()
{
  // PETSc Mat object for stiffness matrix needs to be assembled for this
 
  LOG(DEBUG)<<"applyBoundaryConditions for Laplace";
 
  int nDegreesOfFreedom = this->data_.nDegreesOfFreedom();
  
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

template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::
mesh()
{
  return data_.mesh();
}
  
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::
initialize()
{
  LOG(TRACE) << "FiniteElementMethodBase::initialize";
  
  setStiffnessMatrix();
  setRightHandSide();
  data_.finalAssembly();
  
  applyBoundaryConditions();
}
  
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::
run()
{
  initialize();
  solve();
  data_.print();
  
  outputWriterManager_.writeOutput(data_);
}

template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::
solve()
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