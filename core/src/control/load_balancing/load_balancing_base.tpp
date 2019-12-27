#include "control/load_balancing/load_balancing.h"

#include <omp.h>
#include <sstream>

namespace Control
{

template<class TimeStepping>
LoadBalancingBase<TimeStepping>::
LoadBalancingBase(DihuContext context) :
  Runnable(), ::TimeSteppingScheme::TimeSteppingScheme(context["LoadBalancing"]),
  timeSteppingScheme_(this->context_)
{
}

template<class TimeStepping>
void LoadBalancingBase<TimeStepping>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "LoadBalancingBase::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_; timeStepNo++)
  {
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "LoadBalancingBase, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // set timespan for timeSteppingScheme_
    this->timeSteppingScheme_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance the simulation by the specified time span
    timeSteppingScheme_.advanceTimeSpan();

    // check if the dofs can be rebalanced
    rebalance();
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<class TimeStepping>
void LoadBalancingBase<TimeStepping>::
initialize()
{
  LOG(TRACE) << "LoadBalancingBase::initialize()";

  TimeSteppingScheme::TimeSteppingScheme::initialize();
  timeSteppingScheme_.initialize();
}

template<class TimeStepping>
void LoadBalancingBase<TimeStepping>::
run()
{
  initialize();

  advanceTimeSpan();
}

template<class TimeStepping>
void LoadBalancingBase<TimeStepping>::
reset()
{
  timeSteppingScheme_.reset();
}

template<class TimeStepping>
typename LoadBalancingBase<TimeStepping>::Data &LoadBalancingBase<TimeStepping>::
data()
{
  return timeSteppingScheme_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<class TimeStepping>
typename LoadBalancingBase<TimeStepping>::OutputConnectorDataType &LoadBalancingBase<TimeStepping>::
getOutputConnectorData()
{
  return timeSteppingScheme_.getOutputConnectorData();
}

} // namespace
