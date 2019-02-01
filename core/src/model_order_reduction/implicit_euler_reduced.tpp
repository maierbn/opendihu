#include "model_order_reduction/implicit_euler_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"


namespace ModelOrderReduction
{
template<typename TimeSteppingImplicitType>
ImplicitEulerReduced<TimeSteppingImplicitType>::
ImplicitEulerReduced(DihuContext context):
TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>(context,"ImplicitEulerReduced")
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
  double timeSpan = this->endTime_ - this->startTime_;
  
  LOG(DEBUG) << "ReducedOrderImplicitEuler::advanceTimeSpan(), timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
  << " n steps: " << this->numberTimeSteps_;
  
  // debugging output of matrices
  //this->timestepping_.data_->print();
  
  Vec &solution = this->fullTimestepping_.data().solution()->getValuesContiguous();   // vector of all components in struct-of-array order, as needed by CellML
  Vec &redSolution= this->data().solution()->getValuesContiguous();
  
  Mat &basis = this->dataMOR_->basis()->valuesGlobal();
  
  // loop over time steps
  double currentTime = this->startTime_;
  
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0 )
    {
      LOG(INFO) << "Implicit Euler, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;   
    
    VLOG(1) << "starting from solution: " << *this->fullTimestepping_.data().solution();                   
    
    // adjust the full-order rhs vector such that boundary conditions are satisfied
    this->fullTimestepping_.dirichletBoundaryConditions()->applyInRightHandSide(this->fullTimestepping_.data().solution(), this->fullTimestepping_.dataImplicit().boundaryConditionsRightHandSideSummand());
    
    // advance computed value
    // solve A_R*z^{t+1} = z^{t} for z^{t+1} where A_R is the reduced system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem(redSolution, redSolution);
    
    /*
     PetscInt vec_sz;
     PetscScalar val;
     
     VecGetSize(redSolution,&vec_sz);
     std::cout << "vec_sz in ReducedOrderImplicitEuler: " << vec_sz;
     
     for(int i=0; i< vec_sz; i++)
     {
      VecGetValues(redSolution,1,&i,&val);
      LOG(DEBUG) << "redSolution[" << i << "]: " << val;
     }
    */
    
    // transfer to full-order space
    this->MatMultFull(basis,redSolution,solution);       
    
    VLOG(1) << "solution after integration: " << *this->fullTimestepping_.data().solution();
    
    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);
    
    // write the current output values of the full timestepping
    this->fullTimestepping_.outputWriterManager().writeOutput(this->fullTimestepping_.data(), timeStepNo, currentTime);
    
    // write the current output values of the (reduced) timestepping
    this->outputWriterManager().writeOutput(*this->data_, timeStepNo, currentTime);
    
    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    
    this->fullTimestepping_.data().solution()->restoreValuesContiguous();
  }
  
  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename TimeSteppingImplicitType>
void ImplicitEulerReduced<TimeSteppingImplicitType>::
run()
{
  TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::run();
  
}
  
} //namespace