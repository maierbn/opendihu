#include "model_order_reduction/time_stepping_reduced_explicit.h"

#include <Python.h>
#include<petscmat.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"


namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  ExplicitEulerReduced<TimeSteppingExplicitType>::
  ExplicitEulerReduced(DihuContext context):
  TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>(context,"ExplicitEulerReduced")
  {
  }
   
  template<typename TimeSteppingExplicitType>
  void ExplicitEulerReduced<TimeSteppingExplicitType>::
  advanceTimeSpan()
  {
    // compute timestep width
    double timeSpan = this->endTime_ - this->startTime_;
    
    LOG(DEBUG) << "ReducedOrderExplicitEuler::advanceTimeSpan(), timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
      << " n steps: " << this->numberTimeSteps_;       

    // debugging output of matrices
    //this->fullTimestepping_.data().print();
    
    Vec &solution = this->fullTimestepping_.data().solution()->getValuesContiguous();   // vector of all components in struct-of-array order, as needed by CellML
    Vec &increment = this->fullTimestepping_.data().increment()->getValuesContiguous();
    Vec &redSolution = this->data().solution()->valuesGlobal();
    Vec &redIncrement = this->data().increment()->valuesGlobal();
        
    Mat &basis = this->dataMOR_->basis()->valuesGlobal();
    Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal();    
    
    // loop over time steps
    double currentTime = this->startTime_;
    for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
    {
      if (timeStepNo % this->fullTimestepping_.timeStepOutputInterval() == 0)
      {
        std::stringstream threadNumberMessage;
        threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
        LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
      }
                 
      // full state recovery
      //required in case of operator splitting because only the reduced solutions is transferred.
      this->MatMultFull(basis, redSolution, solution);
            
      VLOG(1) << "starting from full-order solution: " << *this->fullTimestepping_.data().solution();     
      
      // advance computed value
      // compute next delta_u = f(u)
      this->evaluateTimesteppingRightHandSideExplicit(solution, increment, timeStepNo, currentTime);      
      
      VLOG(2) << "computed full-order increment: " << *this->fullTimestepping_.data().increment() << ", dt=" << this->timeStepWidth_;             
      
      // reduction step
      // solution may has been changed inside evaluateTimesteppingRightHandSideExplicit in case of 
      // the stimulation in electrophysiology examples. Therefore, the reduced solution has to be updated.
      this->MatMultReduced(basisTransp, solution, redSolution);
      
      VLOG(2) << "reduced solution before adding the reduced increment" << *this->data().solution();
      
      // reduction of increment
      // modified version of MatMult for MOR
      this->MatMultReduced(basisTransp, increment, redIncrement);
      
      VLOG(2) << "reduced increment: " << *this->data().increment() << ", dt=" << this->timeStepWidth_;             
      
      // integrate, z += dt * delta_z
      VecAXPY(redSolution, this->timeStepWidth_, redIncrement);
      
      VLOG(2) << "reduced solution after integration" << *this->data().solution();           
      
      // advance simulation time
      timeStepNo++;
      currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
      
      // write the current output values of the full-order timestepping
      // full state recovery

      this->MatMultFull(basis, redSolution , solution);
      VLOG(1) << "solution after integration" << *this->fullTimestepping_.data().solution(); 
      
      this->fullTimestepping_.outputWriterManager().writeOutput(this->fullTimestepping_.data(), timeStepNo, currentTime);
      
      // write the current output values of the (reduced) timestepping
      this->outputWriterManager().writeOutput(*this->data_, timeStepNo, currentTime);
                  
    }
    
    this->fullTimestepping_.data().solution()->restoreValuesContiguous();
    this->fullTimestepping_.data().increment()->restoreValuesContiguous();
  }
  
  template<typename TimeSteppingExplicitType>
  void ExplicitEulerReduced<TimeSteppingExplicitType>::
  run()
  {    
    TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>::run();
  }
  
} //namespace
