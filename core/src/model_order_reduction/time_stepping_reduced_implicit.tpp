#include "model_order_reduction/time_stepping_reduced_implicit.h"

#include <Python.h>
#include "utility/python_utility.h"

namespace ModelOrderReduction
{
template<typename TimeSteppingImplicitType>
TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
TimeSteppingSchemeOdeReducedImplicit(DihuContext context):
TimeSteppingSchemeOdeReduced<TimesteppingImplicitType>(context,"ImplicitEuler"), initialized_(false)
{
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
initialize()
{
  if (initialized_)
    return;
  
  LOG(TRACE) << "TimeSteppingSchemeOdeReducedImplicit::initialize()";
  
  // TO BE IMPLEMENTED
  
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->timestepping_.endTime_ - this->timestepping_.startTime_;
  
  LOG(DEBUG) << "ReducedOrderImplicitEuler::advanceTimeSpan(), timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timestepping_.timeStepWidth_
  << " n steps: " << this->timestepping_.numberTimeSteps_;
  
  // debugging output of matrices
  //this->timestepping_.data_->print();
  
  Vec &solution = this->timestepping_.data_->solution().getValuesContiguous();   // vector of all components in struct-of-array order, as needed by CellML
  Vec &redSolution= this->solution().getValuesContiguous();
  
  // loop over time steps
  double currentTime = this->timestepping_.startTime_;
  for (int timeStepNo = 0; timeStepNo < this->timestepping_.numberTimeSteps_;)
  {
    if (timeStepNo % this->timestepping_.timeStepOutputInterval_ == 0)
    {
      std::stringstream threadNumberMessage;
      threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
      LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->timestepping_.numberTimeSteps_<< ", t=" << currentTime;
    }
    
    VLOG(1) << "starting from solution: " << this->timestepping_.data_->solution();
    
    ///
    ///   TO IMPLEMENT
    ///
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    VLOG(1) << "solution after integration: " << this->data_->solution();
    
    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);
  }
  
  this->data_->solution().restoreValuesContiguous();
  
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
run()
{
  TimeSteppingSchemeOdeReduced::run();
  
}
  
} //namespace
