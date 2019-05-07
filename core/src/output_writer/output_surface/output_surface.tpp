#include "output_writer/output_surface/output_surface.h"

namespace OutputWriter
{

template<typename Solver>
OutputSurface<Solver>::
OutputSurface(DihuContext context) :
  context_(context["OutputSurface"]), solver_(context_), data_(context_)
{

}

template<typename Solver>
void OutputSurface<Solver>::
advanceTimeSpan()
{
  solver_.advanceTimeSpan();

  LOG(DEBUG) << "OutputSurface: writeOutput";
  outputWriterManager_.writeOutput(data_);
}

template<typename Solver>
void OutputSurface<Solver>::
initialize()
{
  // initialize output writers
  PythonConfig specificSettings = context_.getPythonConfig();
  outputWriterManager_.initialize(context_, specificSettings);

  // initialize solvers
  solver_.initialize();

  data_.setData(solver_.data());
  data_.initialize();
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
  return solver_.run();
}

template<typename Solver>
void OutputSurface<Solver>::
reset()
{
  return solver_.reset();
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
