#include "time_stepping_scheme/heun.h"

#include <Python.h>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
Heun<DiscretizableInTime>::Heun(DihuContext context) : 
  TimeSteppingSchemeOde<DiscretizableInTime>(context, "Heun")
{
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "Heun");
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename DiscretizableInTime>
void Heun<DiscretizableInTime>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  double timeStepWidth = timeSpan / this->numberTimeSteps_;
  
  // additional storage needed for delta u* = f(u*): do this with petsc!
  //std::vector<double> delta_u_star(this->data_->nUnknowns());
  
 
  LOG(DEBUG) << "Heun::advanceTimeSpan, timeSpan="<<timeSpan<<", timeStepWidth="<<timeStepWidth
    <<" n steps: "<<this->numberTimeSteps_;
  
  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
     LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;
    
    //LOG(DEBUG) << "solution before integration: " << PetscUtility::getStringVector(this->data_->solution().values());
    
    // advance solution value to compute u* first
    // compute next delta_u = f(u_{t})
    this->discretizableInTime_.evaluateTimesteppingRightHandSide(
      this->data_->solution().values(), this->data_->increment().values(), timeStepNo, currentTime);
    
    // integrate u* += dt * delta_u : values = solution.values + timeStepWidth * increment.values 
    VecAXPY(this->data_->solution().values(), timeStepWidth, this->data_->increment().values());
    
    // now, advance solution value to compute u_{t+1}
    //this->discretizableInTime_.evaluateTimesteppingRightHandSide(
    //  this->data_->solution().values(), delta_u_star, timeStepNo + 1, currentTime + timeStepWidth);                  // @ Benni: ist das richtig mit "timeStepNo + 1" und "currentTime + timeStepWidth"?
    
    // integrate u_{t+1} = u_{t} + dt*0.5(delta_u + delta_u_star)
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    //LOG(DEBUG) << "solution after integration: " << PetscUtility::getStringVector(this->data_->solution().values());
    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);
    
    //this->data_->print();
  }
}

template<typename DiscretizableInTime>
void Heun<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}
} // namespace TimeSteppingScheme