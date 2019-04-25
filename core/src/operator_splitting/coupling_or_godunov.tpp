#include "operator_splitting/coupling_or_godunov.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping/time_stepping.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
CouplingOrGodunov<TimeStepping1,TimeStepping2>::
CouplingOrGodunov(DihuContext context, std::string name) :
  OperatorSplitting<TimeStepping1,TimeStepping2>(context, name)  // name is "GodunovSplitting" for Godunov and "Coupling" for Coupling
{
}

template<typename TimeStepping1, typename TimeStepping2>
void CouplingOrGodunov<TimeStepping1,TimeStepping2>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "  CouplingOrGodunov::advanceTimeSpan: timeSpan=[" << this->startTime_<< "," << this->endTime_<< "]"
    << ", n steps: " << this->numberTimeSteps_<< ", timeStepWidth=" << this->timeStepWidth_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << this->schemeName_ << ", timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    LOG(DEBUG) << "  CouplingOrGodunov: time step " << timeStepNo << ", t: " << currentTime;

    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping1 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping1 advanceTimeSpan";

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan();
    
    LOG(DEBUG) << "  CouplingOrGodunov: transfer timeStepping1 -> timeStepping2";
    typename TimeStepping1::TransferableSolutionDataType solutionTimeStepping1 = this->timeStepping1_.getSolutionForTransfer();

    if (VLOG_IS_ON(1))
      VLOG(1) << "  timeStepping1_.getSolutionForTransfer(): " << this->timeStepping1_.getString(solutionTimeStepping1);
    // scale solution in timeStepping1 and transfer to timestepping2_
    SolutionVectorMapping<typename TimeStepping1::TransferableSolutionDataType, typename TimeStepping2::TransferableSolutionDataType>::
      transfer(solutionTimeStepping1, this->timeStepping2_.getSolutionForTransfer());

    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping2 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping2
    this->timeStepping2_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping2 advanceTimeSpan";
    // advance simulation by time span
    this->timeStepping2_.advanceTimeSpan();

    LOG(DEBUG) << "  CouplingOrGodunov: transfer timeStepping2 -> timeStepping1";
    // scale solution in timeStepping2 and transfer to timestepping1_
    SolutionVectorMapping<typename TimeStepping2::TransferableSolutionDataType, typename TimeStepping1::TransferableSolutionDataType>::
      transfer(this->timeStepping2_.getSolutionForTransfer(), this->timeStepping1_.getSolutionForTransfer());

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // store the current simulation in case the program gets interrupted, then the last time gets logged
    Control::PerformanceMeasurement::setParameter("currentSimulationTime", std::to_string(currentTime));
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}
}  // namespace
