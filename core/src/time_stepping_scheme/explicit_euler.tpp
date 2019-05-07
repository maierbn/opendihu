#include "time_stepping_scheme/explicit_euler.h"

#include <Python.h>  // has to be the first included header
#include <omp.h>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
//#include "/usr/local/home/kraemer/opendihu/dependencies/petsc/src/petsc-master/src/vec/vec/impls/seq/seqcuda/cudavecimpl.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ExplicitEuler<DiscretizableInTime>::ExplicitEuler(DihuContext context) :
  TimeSteppingExplicit<DiscretizableInTime>(context, "ExplicitEuler")
{
  LOG(TRACE) << "Entering ExplicitEuler constructor";
  this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->context_); // create data object for explicit euler
  LOG(TRACE) << "Leaving ExplicitEuler constructor";
  
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "ExplicitEuler::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // debugging output of matrices
  //this->data_->print();

  // get vectors of all components in struct-of-array order, as needed by CellML (i.e. one long vector with [state0 state0 state0 ... state1 state1...]
  Vec &solution = this->data_->solution()->getValuesContiguous();
  Vec &increment = this->data_->increment()->getValuesContiguous();

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    PetscErrorCode ierr{0};
    
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Explicit Euler, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    VLOG(1) << "starting from solution: " << this->data_->solution();

    //LOG(WARNING) << "Ready for attachement in explicit_euler.tpp:53";
    //while (timeSpan==0.001)
    //{sleep(5);}
    
    // advance computed value
    // compute next delta_u = f(u)
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
      solution, increment, timeStepNo, currentTime);

    VLOG(1) << "increment: " << this->data_->increment() << ", dt: " << this->timeStepWidth_;

    // integrate, y += dt * delta_u
    // unfortunately, this can only be used with petsc master branch, GNU compiler and "--with-cuda" (approximately)
    //if (this->specificSettings_.hasKey("gpu"))
    //{
    //  ierr = VecAXPY_SeqCUDA(solution, this->timeStepWidth_, increment); CHKERRV(ierr);
    //}
    //else
    //{
    ierr = VecAXPY(solution, this->timeStepWidth_, increment); CHKERRV(ierr);
    //}

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    VLOG(1) << "solution after integration: " << this->data_->solution();

    // apply the prescribed boundary condition values
    this->applyBoundaryConditions();

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  this->data_->solution()->restoreValuesContiguous();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::run()
{
  LOG(TRACE) << "ExplicitEuler::run start";
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
  LOG(TRACE) << "ExplicitEuler::run end";
}
} // namespace TimeSteppingScheme
