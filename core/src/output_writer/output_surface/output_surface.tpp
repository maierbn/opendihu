#include "output_writer/output_surface/output_surface.h"

namespace OutputWriter
{

template<typename Solver>
OutputSurface<Solver>::
OutputSurface(DihuContext context) :
  context_(context["OutputSurface"]), solver_(context_), data_(context_), ownRankInvolvedInOutput_(true)
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
    outputWriterManager_.writeOutput(data_);
  }
}

template<typename Solver>
void OutputSurface<Solver>::
initialize()
{
  if (initialized_)
    return;

  // initialize output writers
  PythonConfig specificSettings = context_.getPythonConfig();
  LOG(DEBUG) << "OutputSurface: initialize output writers";
  outputWriterManager_.initialize(context_, specificSettings);

  // initialize solvers
  solver_.initialize();

  data_.setData(solver_.data());
  data_.initialize();
  ownRankInvolvedInOutput_ = data_.ownRankInvolvedInOutput();
  initialized_ = true;
}

template<typename Solver>
bool OutputSurface<Solver>::
knowsMeshType()
{
  return solver_.knowsMeshType();
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
typename OutputSurface<Solver>::Data3D &OutputSurface<Solver>::
data()
{
  return solver_.data();
}

template<typename Solver>
typename OutputSurface<Solver>::TransferableSolutionDataType OutputSurface<Solver>::
getSolutionForTransfer()
{
  return solver_.getSolutionForTransfer();
}

}  // namespace OutputWriter
