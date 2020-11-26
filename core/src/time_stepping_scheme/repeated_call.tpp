#include "time_stepping_scheme/repeated_call.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace TimeSteppingScheme
{

template<typename Solver>
RepeatedCall<Solver>::RepeatedCall(DihuContext context) :
  TimeSteppingScheme(context["RepeatedCall"]), solver_(context_)
{
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename Solver>
void RepeatedCall<Solver>::
initialize()
{
  if (initialized_)
    return;
 
  TimeSteppingScheme::initialize();
  LOG(TRACE) << "RepeatedCall::initialize";


  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("RepeatedCall", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // initialize underlying Solver object, also with time step width
  solver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();
}

template<typename Solver>
void RepeatedCall<Solver>::advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "RepeatedCall::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "RepeatedCall, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // set sub time step in solver
    this->solver_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance solver
    this->solver_.advanceTimeSpan();

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename Solver>
void RepeatedCall<Solver>::
run()
{
  initialize();
  advanceTimeSpan();
}

} // namespace
