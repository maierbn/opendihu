#include "operator_splitting/coupling_or_godunov.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping/time_stepping.h"
#include "output_connector_data_transfer/output_connector_data_transfer.h"

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


    // --------------- time stepping 1, time span = [0,dt] -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping1 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping1 advanceTimeSpan";

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping1AdvanceTimeSpan_);
    }

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan();
    
    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping1AdvanceTimeSpan_);
      Control::PerformanceMeasurement::start(this->logKeyTransfer12_);
    }

    // --------------- data transfer 1->2 -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov: transfer timeStepping1 -> timeStepping2";

    // set transfer direction 1->2
    this->outputConnection_.setTransferDirection(true);

    // transfer actual values
    std::shared_ptr<typename TimeStepping1::OutputConnectorDataType> solutionTimeStepping1 = this->timeStepping1_.getOutputConnectorData();

    if (VLOG_IS_ON(1))
      VLOG(1) << "  timeStepping1_.getOutputConnectorData(): " << solutionTimeStepping1;

    // transfer to timestepping2_
    SolutionVectorMapping<typename TimeStepping1::OutputConnectorDataType, typename TimeStepping2::OutputConnectorDataType>::
      transfer(solutionTimeStepping1, this->timeStepping2_.getOutputConnectorData(), this->outputConnection_);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTransfer12_);
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping2AdvanceTimeSpan_);
    }

    // --------------- time stepping 2, time span = [0,dt] -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping2 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping2
    this->timeStepping2_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  CouplingOrGodunov: timeStepping2 advanceTimeSpan";
    // advance simulation by time span
    this->timeStepping2_.advanceTimeSpan();

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping2AdvanceTimeSpan_);
      Control::PerformanceMeasurement::start(this->logKeyTransfer21_);
    }

    // --------------- data transfer 2->1 -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov: transfer timeStepping2 -> timeStepping1";

    // set transfer direction 1->2
    this->outputConnection_.setTransferDirection(false);

    // transfer to timestepping1_
    SolutionVectorMapping<typename TimeStepping2::OutputConnectorDataType, typename TimeStepping1::OutputConnectorDataType>::
      transfer(this->timeStepping2_.getOutputConnectorData(), this->timeStepping1_.getOutputConnectorData(), this->outputConnection_);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTransfer21_);
    }

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
