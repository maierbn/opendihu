#include "operator_splitting/coupling_or_godunov.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping/time_stepping.h"
#include "slot_connection/slot_connector_data_transfer.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
void CouplingOrGodunov<TimeStepping1,TimeStepping2>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\")::advanceTimeSpan: timeSpan=[" << this->startTime_<< "," << this->endTime_<< "]"
    << ", n steps: " << this->numberTimeSteps_<< ", timeStepWidth=" << this->timeStepWidth_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << this->schemeName_ << ", timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): time step " << timeStepNo << ", t: " << currentTime;


    // --------------- time stepping 1, time span = [0,dt] -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): timeStepping1 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): timeStepping1 advanceTimeSpan";

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping1AdvanceTimeSpan_);
    }

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan(withOutputWritersEnabled);
    
    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping1AdvanceTimeSpan_);
      Control::PerformanceMeasurement::start(this->logKeyTransfer12_);
    }

    // --------------- data transfer 1->2 -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): transfer timeStepping1 -> timeStepping2";

    // set transfer direction 1->2
    this->slotsConnection_->setTransferDirection(true);

    if (VLOG_IS_ON(1))
      VLOG(1) << "  before transfer 1->2 timeStepping1_.getSlotConnectorData(): " << this->timeStepping1_.getSlotConnectorData();

    // transfer to timestepping2_
    SlotConnectorDataTransfer<typename TimeStepping1::SlotConnectorDataType, typename TimeStepping2::SlotConnectorDataType>::
      transfer(this->timeStepping1_.getSlotConnectorData(), this->timeStepping2_.getSlotConnectorData(), *this->slotsConnection_);

    if (VLOG_IS_ON(1))
      VLOG(1) << "  after transfer 1->2 timeStepping2_.getSlotConnectorData(): " << this->timeStepping2_.getSlotConnectorData();

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTransfer12_);
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping2AdvanceTimeSpan_);
    }

    // --------------- time stepping 2, time span = [0,dt] -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): timeStepping2 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping2
    this->timeStepping2_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): timeStepping2 advanceTimeSpan";
    // advance simulation by time span
    this->timeStepping2_.advanceTimeSpan(withOutputWritersEnabled);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping2AdvanceTimeSpan_);
      Control::PerformanceMeasurement::start(this->logKeyTransfer21_);
    }

    // --------------- data transfer 2->1 -------------------------
    LOG(DEBUG) << "  CouplingOrGodunov(\"" << this->description_ << "\"): transfer timeStepping2 -> timeStepping1";

    // set transfer direction 2->1
    this->slotsConnection_->setTransferDirection(false);

    if (VLOG_IS_ON(1))
      VLOG(1) << "  before transfer 2->1: timeStepping2_.getSlotConnectorData(): " << this->timeStepping2_.getSlotConnectorData();

    // transfer to timestepping1_
    SlotConnectorDataTransfer<typename TimeStepping2::SlotConnectorDataType, typename TimeStepping1::SlotConnectorDataType>::
      transfer(this->timeStepping2_.getSlotConnectorData(), this->timeStepping1_.getSlotConnectorData(), *this->slotsConnection_);

    if (VLOG_IS_ON(1))
      VLOG(1) << "  after transfer 2->1: timeStepping1_.getSlotConnectorData(): " << this->timeStepping1_.getSlotConnectorData();


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
