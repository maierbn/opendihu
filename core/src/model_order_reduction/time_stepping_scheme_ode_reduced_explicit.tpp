#pragma once 

#include "model_order_reduction/time_stepping_scheme_ode_reduced_explicit.h"

#include <Python.h>
#include<petscmat.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"


namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  TimeSteppingSchemeOdeReducedExplicit(DihuContext context):
  TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>(context), initialized_(false) 
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
  evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
  {
    this->timestepping_.discretizableInTime().evaluateTimesteppingRightHandSideExplicit(input, output, timeStepNo, currentTime);   
  }
   
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  advanceTimeSpan()
  {
    // compute timestep width
    double timeSpan = this->endTime_ - this->startTime_;
    
    LOG(DEBUG) << "ReducedOrderExplicitEuler::advanceTimeSpan(), timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;
    
    // debugging output of matrices
    //this->timestepping_.data().print();
    
    Vec &solution = this->timestepping_.data().solution().getContiguousValuesGlobal();   // vector of all components in struct-of-array order, as needed by CellML
    Vec &increment = this->timestepping_.data().increment().getContiguousValuesGlobal();
    
    // loop over time steps
    double currentTime = this->startTime_;
    for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
    {
      if (timeStepNo % this->timestepping_.timeStepOutputInterval() == 0)
      {
        std::stringstream threadNumberMessage;
        threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
        LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
      }
      
      VLOG(1) << "starting from solution: " << this->timestepping_.data().solution();
      
      // advance computed value
      // compute next delta_u = f(u)
      this->evaluateTimesteppingRightHandSideExplicit(solution, increment, timeStepNo, currentTime);
      
      VLOG(1) << "computed increment: " << this->timestepping_.data().increment() << ", dt=" << this->timeStepWidth_;
      
      PetscErrorCode ierr;
      
      // reduction step
      ierr=MatMult(this->data_->basisTransp(),increment,this->data_->redIncrement()); CHKERRV(ierr); 
      
      // integrate, z += dt * delta_z
      VecAXPY(this->data_->redSolution(), this->timeStepWidth_, this->data_->redIncrement());
      
      // full state recovery
      ierr=MatMult(this->data_->basis(),this->data_->redSolution(),solution); CHKERRV(ierr);
      
      // advance simulation time
      timeStepNo++;
      currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
      
      VLOG(1) << "solution after integration: " << this->timestepping_.data().solution();
      
      // write current output values
      this->outputWriterManager_.writeOutput(this->timestepping_.data(), timeStepNo, currentTime);
    }
    
    this->timestepping_.data().solution().restoreContiguousValuesGlobal();
    
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