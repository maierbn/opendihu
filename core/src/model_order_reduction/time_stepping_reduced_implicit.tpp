#include "model_order_reduction/time_stepping_reduced_implicit.h"

#include <Python.h>
#include "utility/python_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "data_management/time_stepping/time_stepping_implicit.h"

namespace ModelOrderReduction
{
template<typename TimeSteppingImplicitType>
TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
TimeSteppingSchemeOdeReducedImplicit(DihuContext context):
TimeSteppingSchemeOdeReduced<TimesteppingImplicitType>(context,"ImplicitEulerReduced"),
initialized_(false)
{
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
initialize()
{
  if (initialized_)
    return;
  
  LOG(TRACE) << "TimeSteppingSchemeOdeReducedImplicit::initialize()";
  
  // TO BE IMPLEMENTED
  TimeSteppingSchemeOdeReduced::initialize();
  
  // compute the system matrix
  this->setSystemMatrix(this->timeStepWidth_);
  
  // initialize the linear solver that is used for solving the implicit system
  initializeLinearSolver();
  
  // set matrix used for linear system and preconditioner to ksp context
  Mat &redSysMatrix = this->dataMOR_->redSysMatrix()->valuesGlobal();
  assert(this->ksp_);
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*ksp_, redSysMatrix, redSysMatrix); CHKERRV(ierr);
  
  // TO BE IMPLEMENTED
  
  TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::initialize();
  
  initialized_=true;
  
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &redSysMatrix = this->dataMOR_->redSysMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  PetscUtility::checkDimensionsMatrixVector(redSysMatrix, input);
  
  // solve the system, KSPSolve(ksp,b,x)
  ierr = KSPSolve(*ksp_, input, output); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp_, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp_, &residualNorm); CHKERRV(ierr);
  
  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp_, &convergedReason); CHKERRV(ierr);
  
  VLOG(1) << "Linear system of implicit time stepping solved in " << numberOfIterations << " iterations, residual norm " << residualNorm
  << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
initializeLinearSolver()
{ 
  if (linearSolver_ == nullptr)
  {
    LOG(DEBUG) << "Implicit time stepping: initialize linearSolver";
    
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->data_->functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  }
  else 
  {
    VLOG(2) << ": linearSolver_ already set";
  }
}


  
} //namespace
