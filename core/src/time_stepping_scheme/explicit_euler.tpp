#include "time_stepping_scheme/explicit_euler.h"

#include <Python.h>

#include "control/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ExplicitEuler<DiscretizableInTime>::ExplicitEuler(const DihuContext &context) : 
  TimeSteppingSchemeOde<DiscretizableInTime>(context)
{
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::extractDict(topLevelSettings, "ExplicitEuler");
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  double timeStepWidth = timeSpan / this->numberTimeSteps_;
 
  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputFrequency_ == 0)
     LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;
    
    // advance computed value
    // compute next delta_u = f(u)
    this->discretizableInTime.evaluateTimesteppingRightHandSide(
      this->data_.solution(), this->data_.increment(), timeStepNo, currentTime);
    
    // integrate, y += dt * delta_u
    VecAXPY(this->data_.solution(), timeStepWidth, this->data_.increment());
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    // write current output values
    this->context_.writeOutput(this->data_, timeStepNo, currentTime);
    
    this->data_.print();
  }
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}
} // namespace TimeSteppingScheme