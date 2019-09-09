#include "output_writer/output_surface/output_surface.h"

namespace OutputWriter
{

template<typename Solver>
OutputSurface<Solver>::
OutputSurface(DihuContext context) :
  context_(context["OutputSurface"]), solver_(context_), data_(context_), ownRankInvolvedInOutput_(true), timeStepNo_(0), currentTime_(0.0)
{

}

template<typename Solver>
void OutputSurface<Solver>::
advanceTimeSpan()
{
  solver_.advanceTimeSpan();

  LOG(DEBUG) << "OutputSurface: writeOutput";
  if (ownRankInvolvedInOutput_)
  {
    outputWriterManager_.writeOutput(data_, timeStepNo_++);
  }
}

template<typename Solver>
void OutputSurface<Solver>::
initialize()
{
  if (initialized_)
    return;


  // initialize solvers
  solver_.initialize();

  data_.setData(solver_.data());
  data_.initialize();
  ownRankInvolvedInOutput_ = data_.ownRankInvolvedInOutput();

  // initialize output writers
  PythonConfig specificSettings = context_.getPythonConfig();
  LOG(DEBUG) << "OutputSurface: initialize output writers";

  // initialize output writer to use smaller rank subset that only contains the ranks that have parts of the surface
  // if the last argument is not given, by default the common rank subset would be used
  if (ownRankInvolvedInOutput_)
  {
    outputWriterManager_.initialize(context_, specificSettings, data_.functionSpace()->meshPartition()->rankSubset());
  }
  initialized_ = true;
}

template<typename Solver>
void OutputSurface<Solver>::
run()
{
  initialize();

  solver_.run();

  LOG(DEBUG) << "OutputSurface: writeOutput";
  if (ownRankInvolvedInOutput_)
  {
    outputWriterManager_.writeOutput(data_);
  }
}

template<typename Solver>
void OutputSurface<Solver>::
reset()
{
  solver_.reset();
}

template<typename Solver>
void OutputSurface<Solver>::
setTimeSpan(double startTime, double endTime)
{
  currentTime_ = startTime;
  solver_.setTimeSpan(startTime, endTime);
}

template<typename Solver>
typename OutputSurface<Solver>::Data3D &OutputSurface<Solver>::
data()
{
  return solver_.data();
}

template<typename Solver>
typename OutputSurface<Solver>::OutputConnectorDataType OutputSurface<Solver>::
getOutputConnectorData()
{
  return solver_.getOutputConnectorData();
}

}  // namespace OutputWriter
