#include "spatial_discretization/finite_element_method/06_timestepping_implicit.h"

#include <Python.h>
#include <iostream>
#include <petscmat.h>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"


namespace SpatialDiscretization
{

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
FiniteElementMethodTimeSteppingImplicit(DihuContext context)
  : FiniteElementMethodTimeSteppingExplicit<BasisOnMeshType, QuadratureType, Term>(context)
{
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
initialize()
{
  this->data_.initialize();
  this->setStiffnessMatrix();
  this->setMassMatrix();
  this->data_.finalAssembly();

  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "You have specified output writers for a FiniteElementMethod which is used for a time stepping problem. "
      "The output will not contain any solution data, only geometry. Probably you want to get output from the time stepping scheme, then define the output writers there.";
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
constexpr int FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
nComponents()
{
  return 1;   // this may be different for structural mechanics
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
Mat &FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
systemMatrix()
{
  return this->systemMatrix_;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
preComputeSystemMatrix(double timeStepWidth)
{
  Mat &systemMatrix=this->systemMatrix_;
  Mat &massMatrix=this->data_.massMatrix;
  Mat &stiffnessMatrix=this->data_.stiffnessMatrix;
  
  
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
solveLinearSystem(Vec &input, Vec &output)
{
  //Au^(t+1)=u^(t)
  
  Mat &systemMatrix=this.systemMatrix;
  
  PetscErrorCode ierr;
  
  // create linear solver context
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp, systemMatrix, systemMatrix); CHKERRV(ierr);
  
   // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(output, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
  
  // solve the system
  ierr = KSPSolve(*ksp, input, output); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);

  //LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  VLOG(1) << "Solution of the linear system recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  
}

/*
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
recoverRightHandSideStrongForm(Vec &result)
{
  // massMatrix * f_strong = rhs_weak
  Vec &rhs = this->data_.rightHandSide().values();   // rhs in weak formulation
  Mat &massMatrix = this->data_.massMatrix();

  PetscErrorCode ierr;

  // create linear solver context
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();

  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp, massMatrix, massMatrix); CHKERRV(ierr);

  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(result, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);

  // solve the system
  ierr = KSPSolve(*ksp, rhs, result); CHKERRV(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);

  //LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  VLOG(1) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
}
*/



template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
checkDimensions(Mat &systemMatrix, Vec &input)
{
#ifndef NDEBUG
  int nRows, nColumns;
  MatGetSize(systemMatrix, &nRows, &nColumns);
  int nEntries;
  VecGetSize(input, &nEntries);
  if (nColumns != nEntries)
  {
    LOG(ERROR) << "Stiffness matrix dimension " << nRows << "x" << nColumns << " does not match input vector (" << nEntries << ")!";
  }
  assert(nColumns == nEntries);
#endif
}

/*
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  // this method computes output = M^{-1}*K*input
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  Vec &rhs = this->data_.rightHandSide().values();

  // check if matrix and vector sizes match
  checkDimensions(stiffnessMatrix, input);

  // compute rhs = stiffnessMatrix*input
  MatMult(stiffnessMatrix, input, rhs);

  // compute output = massMatrix^{-1}*rhs
  recoverRightHandSideStrongForm(output);

  this->data_.print();
  this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);
}
*/

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
knowsMeshType()
{
  return true;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodTimeSteppingImplicit<BasisOnMeshType, QuadratureType, Term>::
mesh()
{
  return FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::mesh();
}

} // namespace SpatialDiscretization