#include "model_order_reduction/time_stepping_scheme_ode_reduced_explicit.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "time_stepping_scheme/time_stepping_scheme.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  TimeSteppingSchemeOdeReducedExplicit(DihuContext context):
  MORBase(context["ModelOrderReduction"]), TimeSteppingScheme(context["ModelOrderReduction"]),
  timestepping_(context_["ModelOrderReduction"]), initialized_(false) 
  {
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  initialize()
  {
    if (initialized_)
      return;
    
    LOG(TRACE) << "TimeSteppingSchemeOdeReducedExplicit::initialize()";
    
    // TO IMPLEMENT 
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  advanceTimeSpan()
  {
    // compute timestep width
    double timeSpan = this->timestepping_.endTime_ - this->timestepping_.startTime_;
    
    LOG(DEBUG) << "ReducedOrderExplicitEuler::advanceTimeSpan(), timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timestepping_.timeStepWidth_
    << " n steps: " << this->timestepping_.numberTimeSteps_;
    
    // debugging output of matrices
    //this->timestepping_.data_->print();
    
    Vec &solution = this->timestepping_.data_->solution().getContiguousValuesGlobal();   // vector of all components in struct-of-array order, as needed by CellML
    Vec &increment = this->timestepping_.data->increment().getContiguousValuesGlobal();
    Vec &redSolution= this->solution().getContiguousValuesGlobal();
    
    // loop over time steps
    double currentTime = this->timestepping_.startTime_;
    for(int timeStepNo = 0; timeStepNo < this->timestepping_.numberTimeSteps_;)
    {
      if (timeStepNo % this->timestepping_.timeStepOutputInterval_ == 0)
      {
        std::stringstream threadNumberMessage;
        threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
        LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->timestepping_.numberTimeSteps_<< ", t=" << currentTime;
      }
      
      VLOG(1) << "starting from solution: " << this->timestepping_.data_->solution();
      
      // advance computed value
      // compute next delta_u = f(u)
      this->timestepping_.discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        solution, increment, timeStepNo, currentTime);
      
      VLOG(1) << "computed increment: " << this->timestepping_.data_->increment() << ", dt=" << this->timestepping_.timeStepWidth_;
      
      PetscErrorCode ierr;
      
      // reduction step
      ierr=MatMult(this->basisTransp(),increment,this->redIncrement_); 
      
      // integrate, y += dt * delta_u
      VecAXPY(this->redSolution_, this->timeStepWidth_, this->redIncrement_);
      
      // full state recovery
      ierr=MatMult(this->basis(),this->redSolution_,solution); 
      
      VLOG(1) << "updated solution: " << this->data_->solution();
      
      // advance simulation time
      timeStepNo++;
      currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
      
      VLOG(1) << "solution after integration: " << this->data_->solution();
      
      // write current output values
      this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);
    }
    
    this->data_->solution().restoreContiguousValuesGlobal();
    
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  run()
  {
    // initialize
    this->initialize();
    
    // do simulations
    this->advanceTimeSpan();
    
  }
  
} //namespace