#include "operator_splitting/godunov.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
Godunov<TimeStepping1,TimeStepping2>::
Godunov(DihuContext context) :
  OperatorSplitting<TimeStepping1,TimeStepping2>(context, "GodunovSplitting")
{
}

template<typename TimeStepping1, typename TimeStepping2>
void Godunov<TimeStepping1,TimeStepping2>::
advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "  Godunov::advanceTimeSpan: timeSpan=[" << this->startTime_<< "," << this->endTime_<< "]"
    << ", n steps: " << this->numberTimeSteps_<< ", timeStepWidth=" << this->timeStepWidth_;

  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Godunov, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    LOG(DEBUG) << "  Godunov: time step " << timeStepNo << ", t: " << currentTime;

    LOG(DEBUG) << "  Godunov: timeStepping1 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping1
    this->timeStepping1_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  Godunov: timeStepping1 advanceTimeSpan";

    // advance simulation by time span
    this->timeStepping1_.advanceTimeSpan();
    
    LOG(DEBUG) << "  Godunov: transfer timeStepping1 -> timeStepping2";
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    this->timeStepping1_.solutionVectorMapping().transfer(*this->timeStepping1_.solution(),
      this->timeStepping2_.solutionVectorMapping(), *this->timeStepping2_.solution());

    LOG(DEBUG) << "  Godunov: timeStepping2 setTimeSpan [" << currentTime << ", " << currentTime+this->timeStepWidth_<< "]";
    // set timespan for timestepping2
    this->timeStepping2_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    LOG(DEBUG) << "  Godunov: timeStepping2 advanceTimeSpan";
    // advance simulation by time span
    this->timeStepping2_.advanceTimeSpan();

    LOG(DEBUG) << "  Godunov: transfer timeStepping2 -> timeStepping1";
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    this->timeStepping2_.solutionVectorMapping().transfer(*this->timeStepping2_.solution(),
      this->timeStepping1_.solutionVectorMapping(), *this->timeStepping1_.solution());

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    LOG(DEBUG) << "  Godunov: write output";
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
