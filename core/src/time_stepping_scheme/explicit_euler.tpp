#include "time_stepping_scheme/explicit_euler.h"

#include <Python.h>  // has to be the first included header
#include <omp.h>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ExplicitEuler<DiscretizableInTime>::ExplicitEuler(DihuContext context) :
  TimeSteppingSchemeOde<DiscretizableInTime>(context, "ExplicitEuler")
{
  this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(context); // create data object for explicit euler
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ExplicitEuler");
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "ExplicitEuler::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // debugging output of matrices
  //this->data_->print();

  Vec &solution = this->data_->solution().getContiguousValuesGlobal();   // vector of all components in struct-of-array order, as needed by CellML
  Vec &increment = this->data_->increment().getContiguousValuesGlobal();

  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
    {
      std::stringstream threadNumberMessage;
      threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
      LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    VLOG(1) << "starting from solution: " << this->data_->solution();

    // advance computed value
    // compute next delta_u = f(u)
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
      solution, increment, timeStepNo, currentTime);

    VLOG(1) << "computed increment: " << this->data_->increment() << ", dt=" << this->timeStepWidth_;

    // integrate, y += dt * delta_u
    VecAXPY(solution, this->timeStepWidth_, increment);

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

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}
} // namespace TimeSteppingScheme
