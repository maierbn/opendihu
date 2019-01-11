#include "model_order_reduction/implicit_euler_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"

namespace ModelOrderReduction
{
template<typename TimeSteppingImplicitType>
ImplicitEulerReduced<TimeSteppingImplicitType>::
ImplicitEulerReduced(DihuContext context):
TimeSteppingSchemeOdeReducedImplicit<TimesteppingImplicitType>(context,"ImplicitEuler"), initialized_(false)
{
}

template<typename TimeSteppingImplicitType>
void ImplicitEulerReduced<TimeSteppingImplicitType>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);
  
  // compute timestep width
  double timeSpan = this->timestepping_.endTime_ - this->timestepping_.startTime_;
  
  LOG(DEBUG) << "ReducedOrderImplicitEuler::advanceTimeSpan(), timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timestepping_.timeStepWidth_
  << " n steps: " << this->timestepping_.numberTimeSteps_;
  
  // debugging output of matrices
  //this->timestepping_.data_->print();
  
  Vec &solution = this->fullTimestepping_.data_->solution().getValuesContiguous();   // vector of all components in struct-of-array order, as needed by CellML
  Vec &redSolution= this->data().solution().getValuesContiguous();
  
  Mat &basis = this->dataMOR_->basis()->valuesGlobal();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal();
  
  // loop over time steps
  double currentTime = this->.startTime_;
  
  for (int timeStepNo = 0; timeStepNo < this->.numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0 )
    {
      LOG(INFO) << "Implicit Euler, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    
    VLOG(1) << "starting from solution: " << this->fullTimestepping_.data().solution();        
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    ///   TO IMPLEMENT
    
    // adjust the full-order rhs vector such that boundary conditions are satisfied
    this->fullTimestepping_->applyInRightHandSide(this->fullTimestepping_->solution(), this->fullTimestepping_->boundaryConditionsRightHandSideSummand());
    
    // transfer rhs to reduced space
    MatMultReduced(basisTransp,solution,redSolution); 
    
    // advance computed value
    // solve A_R*z^{t+1} = z^{t} for z^{t+1} where A_R is the reduced system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem(redSolution, redSolution);
    
    // transfer to full-order space
    MatMultFull(basis,redSolution,solution);
    
    ///   TO IMPLEMENT
    
    VLOG(1) << "solution after integration: " << this->fullTimestepping_.data().solution();
    
    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);
    
    // write current output values
    this->outputWriterManager_.writeOutput(this->fullTimestepping_.data(), timeStepNo, currentTime);
    
    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    
    //this->fullTimestepping_.data().solution()->restoreValuesContiguous();
  }
  
  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename TimeSteppingImplicitType>
void ImplicitEulerReduced<TimeSteppingImplicitType>::
run()
{
  TimeSteppingSchemeOdeReduced::run();
  
}
  
} //namespace