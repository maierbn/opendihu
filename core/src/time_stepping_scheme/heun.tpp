#include "time_stepping_scheme/heun.h"

#include <Python.h>
#include <memory>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
Heun<DiscretizableInTime>::Heun(DihuContext context) :
  TimeSteppingExplicit<DiscretizableInTime>(context, "Heun")
{
}

template<typename DiscretizableInTime>
void Heun<DiscretizableInTime>::initialize()
{
  LOG(TRACE) << "Heun::initialize";

  this->data_ = std::make_shared<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->context_);  // create data object for heun

  // initialize already writes the first output file
  TimeSteppingSchemeOde<DiscretizableInTime>::initialize();
}

template<typename DiscretizableInTime>
void Heun<DiscretizableInTime>::advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "Heun::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // we need to cast the pointer type to the derived class. Otherwise the additional intermediateIncrement()-method of the class TimeSteppingHeun won't be there:
  std::shared_ptr<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>> dataHeun
    = std::static_pointer_cast<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->data_);

  // get vectors of all components in struct-of-array order, as needed by CellML (i.e. one long vector with [state0 state0 state0 ... state1 state1...]
  Vec &solution = this->data_->solution()->getValuesContiguous();
  Vec &increment = this->data_->increment()->getValuesContiguous();
  Vec &intermediateIncrement = dataHeun->intermediateIncrement()->getValuesContiguous();

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "Heun, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    VLOG(1) << "starting from solution (" << this->data_->solution() << "): " << *this->data_->solution();

    // advance solution value to compute u* first
    // compute  delta_u = f(u_{t})
    // we call f(u_{t}) the "increment"
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
      solution, increment, timeStepNo, currentTime);

    // integrate u* += dt * delta_u : values = solution.values + timeStepWidth * increment.values
    VecAXPY(solution, this->timeStepWidth_, increment);

    VLOG(1) << "increment: " << *this->data_->increment() << ", dt: " << this->timeStepWidth_;

    // now, advance solution value to compute u_{t+1}
    // compute  delta_u* = f(u*)
    // we call f(u*) the "intermediateIncrement"
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
      solution, intermediateIncrement, timeStepNo + 1, currentTime + this->timeStepWidth_);

    // we need       u_{t+1} = u_{t} + dt*0.5*(delta_u + delta_u*)
    // however, use: u_{t+1} = u*    + dt*0.5*(delta_u* - delta_u)     (#)
    // where         u*      = u_{t} + dt*delta_u
    //
    // first calculate (f(u*)-f(u_{t})). to save storage we store into f(u*):
    VecAXPY(intermediateIncrement, -1.0, increment);

    // now compute overall step as described above (#)
    VecAXPY(solution, 0.5*this->timeStepWidth_, intermediateIncrement);

    // apply the prescribed boundary condition values
    this->applyBoundaryConditions();

    VLOG(1) << "final solution (" << this->data_->solution() << "): " << *this->data_->solution();

    // check if the solution contains Nans or Inf values
    this->checkForNanInf(timeStepNo, currentTime);

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename DiscretizableInTime>
void Heun<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}

} // namespace TimeSteppingScheme
