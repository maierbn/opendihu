#include "spatial_discretization/finite_element_method/05_time_stepping.h"

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
Mat &FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
systemMatrix()
{
  return this->systemMatrix_;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
preComputeSystemMatrix(double timeStepWidth)
{
  Mat &systemMatrix=this->systemMatrix_;
  Mat &massMatrix=this->data_.massMatrix;
  Mat &stiffnessMatrix=this->data_.stiffnessMatrix;
  
  
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
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

} // namespace SpatialDiscretization