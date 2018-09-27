#include "operator_splitting/strang.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping.h"

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
advanceTimeSpan()
{
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

  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // compute midTime once per step to reuse it. [currentTime, midTime=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    midTime = currentTime + 0.5 * this->timeStepWidth_;

    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Strang, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    LOG(DEBUG) << "  Strang: time step " << timeStepNo << ", t: " << currentTime;

    LOG(DEBUG) << "  Strang: timeStepping1 (first half) setTimeSpan [" << currentTime << ", " << midTime << "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(currentTime, midTime);

    LOG(DEBUG) << "  Strang: timeStepping1 (first half) advanceTimeSpan";

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan();

    LOG(DEBUG) << "  Strang: transfer timeStepping1 -> timeStepping2";
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    this->timeStepping1_.solutionVectorMapping().transfer(*this->timeStepping1_.solution(),
      this->timeStepping2_.solutionVectorMapping(), *this->timeStepping2_.solution());

    LOG(DEBUG) << "  Strang: timeStepping2 (complete) advanceTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping2
    this->timeStepping2_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance simulation by time span
    this->timeStepping2_.advanceTimeSpan();

    LOG(DEBUG) << "  Strang: transfer timeStepping2 -> timeStepping1";
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    this->timeStepping2_.solutionVectorMapping().transfer(*this->timeStepping2_.solution(),
      this->timeStepping1_.solutionVectorMapping(), *this->timeStepping1_.solution());

    LOG(DEBUG) << "  Strang: timeStepping1 (second half) advanceTimeSpan [" << midTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(midTime,currentTime+this->timeStepWidth_);

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan();

    /* option 1. (not implemented)
     * no need to transfer data again, since next operator splitting step will start with timeStepping1 (which has the actual data).
     *
     * Then we could skip output of this->timeStepping2_.data(). However, if they are needed in between two operator splitting steps (for example in the 3-scale muscle model)
     * then this->timeStepping2_.data() should have the actualized data.
     */

    // option 2. transfer data. (implemented)
    LOG(DEBUG) << "  Strang: transfer timeStepping1 -> timeStepping2";
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    this->timeStepping1_.solutionVectorMapping().transfer(*this->timeStepping1_.solution(),
      this->timeStepping2_.solutionVectorMapping(), *this->timeStepping2_.solution());

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    LOG(DEBUG) << "  Strang: write output";
    // write current output values
    if(this->outputData1_)
    {
      LOG(DEBUG) << "write output 1";
      this->outputWriterManager_.writeOutput(this->timeStepping1_.data(), timeStepNo, currentTime);
    }
    if(this->outputData2_)
    {
      LOG(DEBUG) << "write output 2";
      this->outputWriterManager_.writeOutput(this->timeStepping2_.data(), timeStepNo, currentTime);
    }
  }
}

};    // namespace
