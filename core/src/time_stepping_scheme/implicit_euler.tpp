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

  LOG(DEBUG) << "ImplicitEuler::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << timeStepWidth
    << " n steps: " << this->numberTimeSteps_;
  
  int nEntries;
  VecGetSize(this->data_->solution().valuesGlobal(), &nEntries);

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

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // advance computed value
    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix
    this->discretizableInTime_.evaluateTimesteppingRightHandSideImplicit(
      this->data_->solution().valuesGlobal(), this->data_->solution().valuesGlobal(), timeStepNo, currentTime);

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
