#include "time_stepping_scheme/implicit_euler.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
ImplicitEuler<DiscretizableInTimeType>::ImplicitEuler(DihuContext context) :
TimeSteppingImplicit<DiscretizableInTimeType>(context, "ImplicitEuler")
{
}

template<typename DiscretizableInTimeType>
void ImplicitEuler<DiscretizableInTimeType>::advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  
  LOG(DEBUG) << "ImplicitEuler::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  Vec solution = this->data_->solution()->valuesGlobal();

  // loop over time steps
  double currentTime = this->startTime_;
  
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Implicit Euler, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    
    //VLOG(1) << "initial solution: " << *this->data_->solution();

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    // adjust rhs vector such that boundary conditions are satisfied
    this->dirichletBoundaryConditions_->applyInRightHandSide(this->data_->solution(), this->dataImplicit_->boundaryConditionsRightHandSideSummand());

    //VLOG(1) << "solution after apply BC: " << *this->data_->solution();

    // advance computed value
    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem(solution, solution);
    
    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "new solution: " << *this->data_->solution();
    }

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    this->outputWriterManager_.writeOutput(*this->dataImplicit_, timeStepNo, currentTime);
    
    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    //this->data_->print();
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename DiscretizableInTimeType>
void ImplicitEuler<DiscretizableInTimeType>::
setSystemMatrix(double timeStepWidth)
{
  LOG(TRACE) << "setSystemMatrix(timeStepWidth=" << timeStepWidth << ")";

  // compute the system matrix (I - dt*M^{-1}K) where M^{-1} is the lumped mass matrix
  
  Mat &inverseLumpedMassMatrix = this->discretizableInTime_.data().inverseLumpedMassMatrix()->valuesGlobal();
  Mat &stiffnessMatrix = this->discretizableInTime_.data().stiffnessMatrix()->valuesGlobal();
  Mat systemMatrix;
  
  PetscErrorCode ierr;
  
  // compute systemMatrix = M^{-1}K
  // the result matrix is created by MatMatMult
  ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &systemMatrix);
  this->dataImplicit_->initializeSystemMatrix(systemMatrix);
  
  // scale systemMatrix by -dt, systemMatrix = -dt*M^{-1}K
  ierr = MatScale(this->dataImplicit_->systemMatrix()->valuesGlobal(), -timeStepWidth); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I - dt*M^{-1}K
  ierr = MatShift(this->dataImplicit_->systemMatrix()->valuesGlobal(), 1.0); CHKERRV(ierr);
  
  this->dataImplicit_->systemMatrix()->assembly(MAT_FINAL_ASSEMBLY);
  
  //VLOG(1) << *this->dataImplicit_->systemMatrix();
}

} // namespace TimeSteppingScheme
