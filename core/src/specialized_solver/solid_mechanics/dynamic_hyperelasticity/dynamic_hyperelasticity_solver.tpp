#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

namespace TimeSteppingScheme
{

template<typename Term>
DynamicHyperelasticitySolver<Term>::
DynamicHyperelasticitySolver(DihuContext context) :
  TimeSteppingScheme(context["DynamicHyperelasticitySolver"]), staticSolver_(this->context)
{

}

template<typename Term>
void HyperelasticitySolver<Term>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // write reference output values
  this->outputWriterManager_.writeOutput(this->data_, 0, 0.0);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0);

  nonlinearSolve();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
}

template<typename Term>
void HyperelasticitySolver<Term>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

} // namespace TimeSteppingScheme
