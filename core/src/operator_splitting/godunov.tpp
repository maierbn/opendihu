#include "operator_splitting/godunov.h"

#include "control/python_utility.h"
#include "data_management/time_stepping.h"

namespace OperatorSplitting
{
 
template<typename TimeStepping1, typename TimeStepping2>
Godunov<TimeStepping1, TimeStepping2>::
Godunov(const DihuContext &context) :
  OperatorSplitting(context), data_(context),
  timeStepping1_(context_["term1"]), timeStepping2_(context_["term2"])
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::extractDict(topLevelSettings, "OperatorSplitting");
}

template<typename TimeStepping1, typename TimeStepping2>
void Godunov<TimeStepping1, TimeStepping2>::
initialize()
{
  TimeSteppingScheme::initialize();
  timeStepping1_.initialize();
  timeStepping2_.initialize();
  
  outputData1_ = PythonUtility::getOptionBool(specificSettings_, "outputData1", true);
  outputData2_ = PythonUtility::getOptionBool(specificSettings_, "outputData2", true);
}

template<typename TimeStepping1, typename TimeStepping2>
void Godunov<TimeStepping1, TimeStepping2>::
advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  double timeStepWidth = timeSpan / this->numberTimeSteps_;
 
  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_; timeStepNo++)
  {
    LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;
    
    // set timespan for timestepping1
    timeStepping1_.setTimeSpan(currentTime, currentTime+timeStepWidth);
    
    // advance simulation by time span
    timeStepping1_.advanceTimeSpan();
    
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping1_.solutionVectorMapping().transfer(timeStepping1_.solution(), 
      timeStepping2_.solutionVectorMapping(), timeStepping2_.solution());
    
    // set timespan for timestepping2
    timeStepping2_.setTimeSpan(currentTime, currentTime+timeStepWidth);
    
    // advance simulation by time span
    timeStepping2_.advanceTimeSpan();
    
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping2_.solutionVectorMapping().transfer(timeStepping2_.solution(), 
      timeStepping1_.solutionVectorMapping(), timeStepping1_.solution());
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    // write current output values
    if(outputData1_)
      this->context_.writeOutput(timeStepping1_.data(), timeStepNo, currentTime);
    if(outputData2_)
      this->context_.writeOutput(timeStepping2_.data(), timeStepNo, currentTime);
  }
}

template<typename TimeStepping1, typename TimeStepping2>
Vec &Godunov<TimeStepping1, TimeStepping2>::
solution()
{
  return timeStepping1_.solution();
}

template<typename TimeStepping1, typename TimeStepping2>
void Godunov<TimeStepping1, TimeStepping2>::
run()
{
  //1: FE, 2: cellml
  initialize();
  
  advanceTimeSpan();
}

};    // namespace