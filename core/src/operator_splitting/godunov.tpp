#include "operator_splitting/godunov.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping.h"

namespace OperatorSplitting
{
 
template<typename TimeStepping1, typename TimeStepping2>
Godunov<TimeStepping1, TimeStepping2>::
Godunov(const DihuContext &context) :
  OperatorSplitting(context), data_(context),
  timeStepping1_(context_["GodunovSplitting"]["Term1"]),
  timeStepping2_(context_["GodunovSplitting"]["Term2"])
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "GodunovSplitting");
  outputWriterManager_.initialize(specificSettings_);
}

template<typename TimeStepping1, typename TimeStepping2>
void Godunov<TimeStepping1, TimeStepping2>::
initialize()
{
  LOG(TRACE) << "  Godunov::initialize";
  
  TimeSteppingScheme::initialize();
  
  LOG(TRACE) << "  Godunov::initialize done, timeSpan=[" <<this->startTime_<<","<<this->endTime_<<"]"
    <<", n steps: "<<this->numberTimeSteps_;
  
  // initialize time stepping objects, if only one knows its MeshType, initialize that first
  //(e.g. CellML-adapter does not know, because it is independent of the type of a mesh)
  if (timeStepping2_.knowsMeshType() && !timeStepping1_.knowsMeshType())
  {
    LOG(DEBUG) << "  Godunov::initialize timeStepping2";
    timeStepping2_.initialize();
    LOG(DEBUG) << "  Godunov::initialize timeStepping1";
    timeStepping1_.initialize();
  }
  else
  {
    LOG(DEBUG) << "  Godunov::initialize timeStepping1";
    timeStepping1_.initialize();
    LOG(DEBUG) << "  Godunov::initialize timeStepping2";
    timeStepping2_.initialize();
  }
  
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
 
  LOG(DEBUG) << "  Godunov::advanceTimeSpan: timeSpan=[" <<this->startTime_<<","<<this->endTime_<<"]"
    <<", n steps: "<<this->numberTimeSteps_<<", timeStepWidth="<<timeStepWidth;
  
  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;
    LOG(DEBUG) << "  Godunov: time step "<<timeStepNo<<", t: "<<currentTime;
    
    LOG(DEBUG) << "  Godunov: timeStepping1 setTimeSpan ["<<currentTime<<", "<<currentTime+timeStepWidth<<"]";
    // set timespan for timestepping1
    timeStepping1_.setTimeSpan(currentTime, currentTime+timeStepWidth);
    
    LOG(DEBUG) << "  Godunov: timeStepping1 advanceTimeSpan"; 
    
    // advance simulation by time span
    timeStepping1_.advanceTimeSpan();
    
    LOG(DEBUG) << "  Godunov: transfer timeStepping1 -> timeStepping2"; 
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping1_.solutionVectorMapping().transfer(timeStepping1_.solution(), 
      timeStepping2_.solutionVectorMapping(), timeStepping2_.solution());
    
    LOG(DEBUG) << "  Godunov: timeStepping2 setTimeSpan ["<<currentTime<<", "<<currentTime+timeStepWidth<<"]";
    // set timespan for timestepping2
    timeStepping2_.setTimeSpan(currentTime, currentTime+timeStepWidth);
    
    LOG(DEBUG) << "  Godunov: timeStepping2 advanceTimeSpan";
    // advance simulation by time span
    timeStepping2_.advanceTimeSpan();
    
    LOG(DEBUG) << "  Godunov: transfer timeStepping2 -> timeStepping1"; 
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping2_.solutionVectorMapping().transfer(timeStepping2_.solution(), 
      timeStepping1_.solutionVectorMapping(), timeStepping1_.solution());
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    LOG(DEBUG) << "  Godunov: write output"; 
    // write current output values
    if(outputData1_)
    {
      LOG(DEBUG) << "write output 1";
      this->outputWriterManager_.writeOutput(timeStepping1_.data(), timeStepNo, currentTime);
    }
    if(outputData2_)
    {
      LOG(DEBUG) << "write output 2";
      this->outputWriterManager_.writeOutput(timeStepping2_.data(), timeStepNo, currentTime);
    }
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
  // initialize data structurures
  initialize();
  
  // run simulation
  advanceTimeSpan();
}

template<typename TimeStepping1, typename TimeStepping2>
bool Godunov<TimeStepping1, TimeStepping2>::
knowsMeshType()
{
  return timeStepping1_.knowsMeshType() && timeStepping2_.knowsMeshType();
}
  
};    // namespace