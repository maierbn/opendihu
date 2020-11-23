#include "operator_splitting/strang.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping/time_stepping.h"
#include "slot_connection/slot_connector_data.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
Strang<TimeStepping1,TimeStepping2>::
Strang(DihuContext context) :
  OperatorSplitting<TimeStepping1,TimeStepping2>(context, "StrangSplitting")
{
}

template<typename TimeStepping1, typename TimeStepping2>
void Strang<TimeStepping1,TimeStepping2>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "  Strang::advanceTimeSpan: timeSpan=[" << this->startTime_<< "," << this->endTime_<< "]"
    << ", n steps: " << this->numberTimeSteps_<< ", timeStepWidth=" << this->timeStepWidth_;

  // picture for strang splitting:
  // ===== t ==>
  // -1->
  //  /
  // ---2--->
  //      /
  //     -1->
  //        |
  //        2

  // loop over time steps
  double currentTime = this->startTime_;
  double midTime = 0.0;

  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // compute midTime once per step to reuse it. [currentTime, midTime=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    midTime = currentTime + 0.5 * this->timeStepWidth_;

    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "Strang, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    LOG(DEBUG) << "  Strang: time step " << timeStepNo << ", t: " << currentTime;
    LOG(DEBUG) << "  Strang: timeStepping1 (first half) setTimeSpan [" << currentTime << ", " << midTime << "]";

    // --------------- time stepping 1, time span = [0,midTime] -------------------------
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping1AdvanceTimeSpan_);

    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(currentTime, midTime);

    LOG(DEBUG) << "  Strang: timeStepping1 (first half) advanceTimeSpan";

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan(withOutputWritersEnabled);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping1AdvanceTimeSpan_);
      Control::PerformanceMeasurement::start(this->logKeyTransfer12_);
    }

    // --------------- data transfer 1->2 -------------------------
    LOG(DEBUG) << "  Strang: transfer timeStepping1 -> timeStepping2";

    // set transfer direction 1->2
    this->slotsConnection_->setTransferDirection(true);

    // transfer to timestepping2_
    SlotConnectorDataTransfer<typename TimeStepping1::SlotConnectorDataType, typename TimeStepping2::SlotConnectorDataType>::
      transfer(this->timeStepping1_.getSlotConnectorData(), this->timeStepping2_.getSlotConnectorData(), *this->slotsConnection_);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTransfer12_);
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping2AdvanceTimeSpan_);
    }

    // --------------- time stepping 2, time span = [0,dt] -------------------------
    LOG(DEBUG) << "  Strang: timeStepping2 (complete) advanceTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";

    // set timespan for timestepping2
    this->timeStepping2_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance simulation by time span
    this->timeStepping2_.advanceTimeSpan(withOutputWritersEnabled);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping2AdvanceTimeSpan_);
      Control::PerformanceMeasurement::start(this->logKeyTransfer21_);
    }

    // --------------- data transfer 2->1 -------------------------
    LOG(DEBUG) << "  Strang: transfer timeStepping2 -> timeStepping1";

    // set transfer direction 2->1
    this->slotsConnection_->setTransferDirection(false);

    // scale solution in timeStepping2 and transfer to timestepping1_
    SlotConnectorDataTransfer<typename TimeStepping2::SlotConnectorDataType, typename TimeStepping1::SlotConnectorDataType>::
      transfer(this->timeStepping2_.getSlotConnectorData(), this->timeStepping1_.getSlotConnectorData(), *this->slotsConnection_);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTransfer21_);
      Control::PerformanceMeasurement::start(this->logKeyTimeStepping1AdvanceTimeSpan_);
    }

    // --------------- time stepping 1, time span = [midTime,dt] -------------------------
    LOG(DEBUG) << "  Strang: timeStepping1 (second half) advanceTimeSpan [" << midTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(midTime,currentTime+this->timeStepWidth_);

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan(withOutputWritersEnabled);

    if (this->durationLogKey_ != "")
    {
      Control::PerformanceMeasurement::stop(this->logKeyTimeStepping1AdvanceTimeSpan_);
    }

    /* option 1. (implemented)
     * no need to transfer data again, since next operator splitting step will start with timeStepping1 (which has the actual data).
     *
     * Then we could skip output of this->timeStepping2_.data(). However, if they are needed in between two operator splitting steps (for example in the 3-scale muscle model)
     * then this->timeStepping2_.data() should have the actualized data.
     */

    /* option 2. transfer data. (not implemented)
     * the second option to not transfer the data again is needed, because two subsequent transfers in the same direction yield an error with shared field variables.
     * When this will be needed in the future, there have to be actually two transfers timeStepping1->timeStepping2 and timeStepping->timeStepping1 again,
     * to have a matching VecResetArray to the previous VecPlaceArray.
     */
#if 0
    LOG(DEBUG) << "  Strang: transfer timeStepping1 -> timeStepping2";

    // set transfer direction 1->2
    this->slotsConnection_->setTransferDirection(true);

    // transfer to timestepping2_
    SlotConnectorDataTransfer<typename TimeStepping1::SlotConnectorDataType, typename TimeStepping2::SlotConnectorDataType>::
      transfer(this->timeStepping1_.getSlotConnectorData(), this->timeStepping2_.getSlotConnectorData(), this->slotsConnection_);
#endif

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
