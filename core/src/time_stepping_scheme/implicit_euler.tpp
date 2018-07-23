#include "time_stepping_scheme/implicit_euler.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ImplicitEuler<DiscretizableInTime>::ImplicitEuler(DihuContext context) :
  TimeSteppingSchemeOde<DiscretizableInTime>(context, "ImplicitEuler")
{
  this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTime::BasisOnMesh, DiscretizableInTime::nComponents()>>(context); // create data object for implicit euler
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ImplicitEuler");
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename DiscretizableInTime>
void ImplicitEuler<DiscretizableInTime>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  double timeStepWidth = timeSpan / this->numberTimeSteps_;

  LOG(DEBUG) << "ExplicitEuler::advanceTimeSpan, timeSpan="<<timeSpan<<", timeStepWidth="<<timeStepWidth
    <<" n steps: "<<this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
     LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;

    // advance computed value
    // compute next delta_u = f(u)
    this->discretizableInTime_.evaluateTimesteppingRightHandSide(
      this->data_->solution().values(), this->data_->increment().values(), timeStepNo, currentTime);

    // integrate, y += dt * delta_u
    VecAXPY(this->data_->solution().values(), timeStepWidth, this->data_->increment().values());

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
void ImplicitEuler<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}
} // namespace TimeSteppingScheme