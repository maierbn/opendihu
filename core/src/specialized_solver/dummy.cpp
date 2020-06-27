#include "specialized_solver/dummy.h"

#include "control/diagnostic_tool/solver_structure_visualizer.h"

Dummy::
Dummy(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context),
  data_(this->context_)
{
}

void Dummy::
advanceTimeSpan()
{
}

void Dummy::initialize()
{
  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("Dummy", false);   // hasInternalConnectionToFirstNestedSolver=false (the last argument) means output connector data is not shared with the first subsolver

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize();

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());
}

void Dummy::run()
{
  initialize();
  advanceTimeSpan();
}

void Dummy::reset()
{
}

typename Dummy::Data &Dummy::data()
{
  // get a reference to the data object
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
std::shared_ptr<typename Dummy::OutputConnectorDataType> Dummy::
getOutputConnectorData()
{
  return data_.getOutputConnectorData();
}
