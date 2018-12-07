#include "model_order_reduction/time_stepping_reduced_explicit.h"

#include <Python.h>
#include<petscmat.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"


namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  TimeSteppingSchemeOdeReducedExplicit(DihuContext context):
  TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>(context,"ExplicitEulerReduced"), initialized_(false) 
  {
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  initialize()
  {
    if (initialized_)
      return;
    
    LOG(TRACE) << "TimeSteppingSchemeOdeReducedExplicit::initialize()";
    
    TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>::initialize(); 
    
    initialized_=true;
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
  {
    this->fullTimestepping_.discretizableInTime().evaluateTimesteppingRightHandSideExplicit(input, output, timeStepNo, currentTime);
    //ISCreateGeneral(comm,n, idx[],mode,&is);
    //VecGetSubVector(X,is,&Y) should it be called every time step?
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
    //this->fullTimestepping_.data().print();
    
    Vec &solution = this->fullTimestepping_.data().solution()->getValuesContiguous();   // vector of all components in struct-of-array order, as needed by CellML
    Vec &increment = this->fullTimestepping_.data().increment()->getValuesContiguous();
    Vec &redSolution= this->data_->solution()->getValuesContiguous();
    Vec &redIncrement= this->data_->increment()->getValuesContiguous();
    
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
      
      VLOG(1) << "starting from solution: " << this->fullTimestepping_.data().solution();
      
      // advance computed value
      // compute next delta_u = f(u)
      this->evaluateTimesteppingRightHandSideExplicit(solution, increment, timeStepNo, currentTime);
      
      VLOG(1) << "computed increment: " << this->fullTimestepping_.data().increment() << ", dt=" << this->timeStepWidth_;
      
      PetscErrorCode ierr;
      PetscInt increment_size, redIncrement_size, mat_sz_1, mat_sz_2;
      
      VecGetSize(increment,&increment_size);
      VecGetSize(redIncrement,&redIncrement_size);
      MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
      
      LOG(DEBUG) << "increment_size: " << increment_size << "========================";
      LOG(DEBUG) << "redIncrement_size: " << redIncrement_size;
      LOG(DEBUG) << "mat_sz_1: " << mat_sz_1 << "mat_sz_2: " << mat_sz_2;
      // reduction step
      ierr=MatMult(basisTransp, increment, redIncrement); CHKERRV(ierr); 
      
      // integrate, z += dt * delta_z
      VecAXPY(redSolution, this->timeStepWidth_, redIncrement);
      
      // full state recovery
      ierr=MatMult(basis, redSolution , solution); CHKERRV(ierr);
      
      // advance simulation time
      timeStepNo++;
      currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
      
      VLOG(1) << "solution after integration: " << this->fullTimestepping_.data().solution();
      
      // write current output values
      this->outputWriterManager_.writeOutput(this->fullTimestepping_.data(), timeStepNo, currentTime);
    }
    
    this->fullTimestepping_.data().solution()->restoreValuesContiguous();
    
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
