#include "time_stepping_scheme/time_stepping_implicit.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
TimeSteppingImplicit<DiscretizableInTimeType>::TimeSteppingImplicit(DihuContext context, std::string name) :
TimeSteppingSchemeOde<DiscretizableInTimeType>(context, name)
{
  this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>(context); // create data object for implicit euler
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, name);
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTimeType>::run();
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
initialize()
{
  if (this->initialized_)
    return;
  
  TimeSteppingSchemeOde<DiscretizableInTimeType>::initialize();
  
  // compute the system matrix
  this->setSystemMatrix(this->timeStepWidth_);
  Mat &systemMatrix = this->data_->systemMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  
  initializeLinearSolver();
  
  // set matrix used for linear system and preconditioner to ksp context
  assert(this->ksp_);
  ierr = KSPSetOperators (*ksp_, systemMatrix, systemMatrix); CHKERRV(ierr);
  
  this->initialized_ = true;
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &systemMatrix = this->data_->systemMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  PetscUtility::checkDimensionsMatrixVector(systemMatrix, input);
  
  // solve the system
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

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
initializeLinearSolver()
{ 
  if (linearSolver_ == nullptr)
  {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    LOG(DEBUG) << s.str() << ": ImplicitEuler: initialize linearSolver";
    
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->data_->functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  }
  else 
  {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    VLOG(2) << s.str() << ": linearSolver_ already set";
    
  }
}

} // namespace TimeSteppingScheme
