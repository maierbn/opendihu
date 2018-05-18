#include "time_stepping_scheme/heun.h"

#include <Python.h>
#include <memory>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
Heun<DiscretizableInTime>::Heun(DihuContext context) :
  TimeSteppingSchemeOde<DiscretizableInTime>(context, "Heun")
{
  this->data_ = std::make_shared <Data::TimeSteppingHeun<typename DiscretizableInTime::BasisOnMesh, DiscretizableInTime::nComponents()>>(context);  // create data object for heun
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

  LOG(DEBUG) << "Heun::advanceTimeSpan, timeSpan="<<timeSpan<<", timeStepWidth="<<timeStepWidth
    <<" n steps: "<<this->numberTimeSteps_;

  // we need to cast the pointer type to the derived class. Otherwise the additional intermediateIncrement()-method of the class TimeSteppingHeun won't be there:
  std::shared_ptr<Data::TimeSteppingHeun<typename DiscretizableInTime::BasisOnMesh, DiscretizableInTime::nComponents()>> dataHeun
    = std::static_pointer_cast<Data::TimeSteppingHeun<typename DiscretizableInTime::BasisOnMesh, DiscretizableInTime::nComponents()>>(this->data_);

  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
     LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;

    //LOG(DEBUG) << "solution before integration: " << PetscUtility::getStringVector(this->data_->solution().values());

    // advance solution value to compute u* first
    // compute  delta_u = f(u_{t})
    // we call f(u_{t}) the "increment"
    this->discretizableInTime_.evaluateTimesteppingRightHandSide(
      this->data_->solution().values(), this->data_->increment().values(), timeStepNo, currentTime);

    // integrate u* += dt * delta_u : values = solution.values + timeStepWidth * increment.values
    VecAXPY(this->data_->solution().values(), timeStepWidth, this->data_->increment().values());

    // now, advance solution value to compute u_{t+1}
    // compute  delta_u* = f(u*)
    // we call f(u*) the "intermediateIncrement"
    this->discretizableInTime_.evaluateTimesteppingRightHandSide(
      this->data_->solution().values(), dataHeun->intermediateIncrement().values(), timeStepNo + 1, currentTime + timeStepWidth);

    // integrate u_{t+1} = u_{t} + dt*0.5(delta_u + delta_u_star)
    // however, use: u_{t+1} = u* + 0.5*dt*(f(u*)-f(u_{t}))     (#)
    //
    // first calculate (f(u*)-f(u_{t})). to save storage we store into f(u*):
    VecAXPY( dataHeun->intermediateIncrement().values(),-1.0,this->data_->increment().values());
    // now compute overall step as described above (#)
    VecAXPY(this->data_->solution().values(), 0.5*timeStepWidth,  dataHeun->intermediateIncrement().values());

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