#include "specialized_solver/my_new_solver/my_new_timestepping_solver.h"

#include <omp.h>
#include <sstream>

MuscleContractionSolver::
MuscleContractionSolver(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context["MuscleContractionSolver"]),
  dynamicHyperelasticitySolver_(this->context_)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

}

void MuscleContractionSolver::
advanceTimeSpan()
{
  // This method computes some time steps of the simulation by running a for loop over the time steps.
  // The number of steps, timestep width and current time are all set by the parent class, TimeSteppingScheme.

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute time span of this method
  double timeSpan = this->endTime_ - this->startTime_;

  // output for debugging
  LOG(DEBUG) << "MuscleContractionSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // in defined intervals (settings "timeStepOutputInterval") print out the current timestep
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "MuscleContractionSolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    this->dynamicHyperelasticitySolver_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance the simulation by the specified time span
    dynamicHyperelasticitySolver_.advanceTimeSpan();

    // advance simulation time
    timeStepNo++;

    // compute new current simulation time
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values using the output writers
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}


void MuscleContractionSolver::
initialize()
{
  // initialize() will be called before the simulation starts.

  // call initialize of the parent class, this parses the timestepping settings from the settings file
  TimeSteppingScheme::TimeSteppingScheme::initialize();

  // call initialize of the nested timestepping solver
  dynamicHyperelasticitySolver_.initialize();
}


void MuscleContractionSolver::
run()
{
  // The run method should not be changed. It is the method that gets called directly from the example main file.
  // If this solver itself is nested in other solvers or coupling schemes,
  // run() will not be called, but the surrounding solver will call initialize() and advanceTimeSpan().
  initialize();

  advanceTimeSpan();
}


void MuscleContractionSolver::
reset()
{
  dynamicHyperelasticitySolver_.reset();

  // "uninitialize" everything
}

typename MuscleContractionSolver::Data &MuscleContractionSolver::
data()
{
  // get a reference to the data object of dynamicHyperelasticitySolver_
  return dynamicHyperelasticitySolver_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class

std::shared_ptr<typename MuscleContractionSolver::OutputConnectorDataType> MuscleContractionSolver::
getOutputConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the timeSteppingScheme_.
  return dynamicHyperelasticitySolver_.getOutputConnectorData();
}
