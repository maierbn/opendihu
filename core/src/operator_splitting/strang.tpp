#include "operator_splitting/strang.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping.h"

namespace OperatorSplitting
{
 
template<typename TimeStepping1, typename TimeStepping2>
Strang<TimeStepping1, TimeStepping2>::
Strang(const DihuContext &context) :
  OperatorSplitting(context), data_(context),
  timeStepping1_(context_["StrangSplitting"]["Term1"]),
  timeStepping2_(context_["StrangSplitting"]["Term2"])
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "StrangSplitting");
  outputWriterManager_.initialize(specificSettings_);
}

template<typename TimeStepping1, typename TimeStepping2>
void Strang<TimeStepping1, TimeStepping2>::
initialize()
{
  LOG(TRACE) << "  Strang::initialize";
  
  TimeSteppingScheme::initialize();
  
  LOG(TRACE) << "  Strang::initialize done, timeSpan=[" <<this->startTime_<<","<<this->endTime_<<"]"
    <<", n steps: "<<this->numberTimeSteps_;
  
  // initialize time stepping objects, if only one knows its MeshType, initialize that first
  //(e.g. CellML-adapter does not know, because it is independent of the type of a mesh)
  if (timeStepping2_.knowsMeshType() && !timeStepping1_.knowsMeshType())
  {
    LOG(DEBUG) << "  Strang::initialize timeStepping2";
    timeStepping2_.initialize();
    LOG(DEBUG) << "  Strang::initialize timeStepping1";
    timeStepping1_.initialize();
  }
  else
  {
    LOG(DEBUG) << "  Strang::initialize timeStepping1";
    timeStepping1_.initialize();
    LOG(DEBUG) << "  Strang::initialize timeStepping2";
    timeStepping2_.initialize();
  }
  
  outputData1_ = PythonUtility::getOptionBool(specificSettings_, "outputData1", true);
  outputData2_ = PythonUtility::getOptionBool(specificSettings_, "outputData2", true);
}

template<typename TimeStepping1, typename TimeStepping2>
void Strang<TimeStepping1, TimeStepping2>::
advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  double timeStepWidth = timeSpan / this->numberTimeSteps_;
 
  LOG(DEBUG) << "  Strang::advanceTimeSpan: timeSpan=[" <<this->startTime_<<","<<this->endTime_<<"]"
    <<", n steps: "<<this->numberTimeSteps_<<", timeStepWidth="<<timeStepWidth;
  
  // loop over time steps
  double currentTime = this->startTime_;
  double midTime = 0.0;
  
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // compute midTime once per step to reuse it. [currentTime, midTime=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    midTime = currentTime + 0.5 * timeStepWidth;
    
    LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;
    LOG(DEBUG) << "  Strang: time step "<<timeStepNo<<", t: "<<currentTime;
    
    LOG(DEBUG) << "  Strang: timeStepping1 (first half) setTimeSpan ["<<currentTime<<", "<<midTime<<"]";
    // set timespan for timestepping1
    timeStepping1_.setTimeSpan(currentTime, midTime);
    
    LOG(DEBUG) << "  Strang: timeStepping1 (first half) advanceTimeSpan"; 
    
    // advance simulation by time span
    timeStepping1_.advanceTimeSpan();
    
    LOG(DEBUG) << "  Strang: transfer timeStepping1 -> timeStepping2"; 
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping1_.solutionVectorMapping().transfer(timeStepping1_.solution(), 
      timeStepping2_.solutionVectorMapping(), timeStepping2_.solution());
    
    LOG(DEBUG) << "  Strang: timeStepping2 (complete) setTimeSpan ["<<currentTime<<", "<<currentTime+timeStepWidth<<"]";
    // set timespan for timestepping2
    timeStepping2_.setTimeSpan(currentTime, currentTime+timeStepWidth);
    
    LOG(DEBUG) << "  Strang: timeStepping2 (complete) advanceTimeSpan";
    // advance simulation by time span
    timeStepping2_.advanceTimeSpan();
    
    LOG(DEBUG) << "  Strang: transfer timeStepping2 -> timeStepping1"; 
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping2_.solutionVectorMapping().transfer(timeStepping2_.solution(), 
      timeStepping1_.solutionVectorMapping(), timeStepping1_.solution());
    
    LOG(DEBUG) << "  Strang: timeStepping1 (second half) setTimeSpan ["<<midTime<<", "<<currentTime+timeStepWidth<<"]";
    // set timespan for timestepping1
    timeStepping1_.setTimeSpan(midTime,currentTime+timeStepWidth);
    
    LOG(DEBUG) << "  Strang: timeStepping1 (second half) advanceTimeSpan"; 
    
    // advance simulation by time span
    timeStepping1_.advanceTimeSpan();
    
    /* option 1. (not implemented)
     * no need to transfer data again, since next operator splitting step will start with timeStepping1 (which has the actual data).
     * 
     * Then we could skip output of timeStepping2_.data(). However, if they are needed in between two operator splitting steps (for example in the 3-scale muscle model)
     * then timeStepping2_.data() should have the actualized data.
     */
    
    // option 2. transfer data. (implemented)
    LOG(DEBUG) << "  Strang: transfer timeStepping1 -> timeStepping2"; 
    // transfer data from timestepping1_.data_.solution_ to timestepping2_.data_.solution_
    timeStepping1_.solutionVectorMapping().transfer(timeStepping1_.solution(), 
      timeStepping2_.solutionVectorMapping(), timeStepping2_.solution());
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    LOG(DEBUG) << "  Strang: write output"; 
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
Vec &Strang<TimeStepping1, TimeStepping2>::
solution()
{
  return timeStepping1_.solution();
}

template<typename TimeStepping1, typename TimeStepping2>
void Strang<TimeStepping1, TimeStepping2>::
run()
{
  // initialize data structurures
  initialize();
  
  // run simulation
  advanceTimeSpan();
}

template<typename TimeStepping1, typename TimeStepping2>
bool Strang<TimeStepping1, TimeStepping2>::
knowsMeshType()
{
  return timeStepping1_.knowsMeshType() && timeStepping2_.knowsMeshType();
}
  
};    // namespace