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
TimeSteppingImplicit<DiscretizableInTimeType>::TimeSteppingImplicit(DihuContext context) :
TimeSteppingSchemeOde<DiscretizableInTimeType>(context, "ImplicitEuler")
{
  this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>(context); // create data object for implicit euler
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ImplicitEuler");
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "ImplicitEuler::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;
  
  int nEntries;
  VecGetSize(this->data_->solution().valuesGlobal(), &nEntries);

  // loop over time steps
  double currentTime = this->startTime_;
  
  this->setSystemMatrix(this->timeStepWidth_);
  Mat &systemMatrix = this->data_.systemMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  
  // set matrix used for linear system and preconditioner to ksp context
  assert(this->ksp_);
  ierr = KSPSetOperators (*ksp_, systemMatrix, systemMatrix); CHKERRV(ierr);
  
  
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
    {
      std::stringstream threadNumberMessage;
      threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
      LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // advance computed value
    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix
    
    this->solveLinearSystem(this->data_->solution().valuesGlobal(), this->data_->solution().valuesGlobal());

    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

    //this->data_->print();
  }
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTimeType>::run();
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &systemMatrix = this->data_.systemMatrix()->valuesGlobal();
  
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

/*
template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
setSystemMatrix(double timeStepWidth)
{
  LOG(TRACE) << "setSystemMatrix(timeStepWidth=" << timeStepWidth << ")";
  
  if(!this->discretizableInTime_.invLumMassMatrixSet())
    this->discretizableInTime_.setInvLumMassMatrix();
  
  // compute the system matrix (I - dt*M^{-1}K) where M^{-1} is the lumped mass matrix
  
  Mat &inverseLumpedMassMatrix = this->data_.inverseLumpedMassMatrix()->valuesGlobal();
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  
  // compute systemMatrix = M^{-1}K
  // the result matrix is created by MatMatMult
  ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &systemMatrix_);
  this->data_.initializeSystemMatrix(systemMatrix);
  
  // scale systemMatrix by -dt, systemMatrix = -dt*M^{-1}K
  ierr = MatScale(this->data_.systemMatrix()->valuesGlobal(), -timeStepWidth); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I - dt*M^{-1}K
  ierr = MatShift(this->data_.systemMatrix()->valuesGlobal(), 1.0); CHKERRV(ierr);
  
  systemMatrix_->assembly(MAT_FINAL_ASSEMBLY);
  
  VLOG(1) << *this->data_.systemMatrix();
}
*/

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
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
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
