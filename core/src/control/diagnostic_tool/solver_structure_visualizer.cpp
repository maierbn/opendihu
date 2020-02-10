#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "data_management/output_connector_data.h"
#include "output_connector_data_transfer/output_connection.h"

//! constructor
SolverStructureVisualizer::SolverStructureVisualizer()
{
  // initialize root solver
  solverRoot_ = std::make_shared<solver_t>();
  solverRoot_->parent = solverRoot_;

  currentSolver_ = solverRoot_;
  enabled_ = true;
  nDisableCalls_ = 0;
}

void SolverStructureVisualizer::addSolver(std::string name)
{
  LOG(DEBUG) << "SolverStructureVisualizer::addSolver(\"" <<  name << "\") under \"" << currentSolver_->parent->name << "\", nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_  << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  if (!currentSolver_)
    LOG(FATAL) << "addSolver on invalid currentSolver";
  currentSolver_->name = name;

  LOG(DEBUG) << "addSolver \"" << name << "\".";
}

//! indicate that all further calls to addSolver will be children of the current solver
void SolverStructureVisualizer::beginChild(std::string description)
{
  LOG(DEBUG) << "SolverStructureVisualizer::beginChild() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  if (!currentSolver_)
    LOG(FATAL) << "beginChild on invalid currentSolver";

  std::shared_ptr<solver_t> newChild = std::make_shared<solver_t>();
  newChild->parent = currentSolver_;
  newChild->description = description;
  currentSolver_->children.push_back(newChild);

  LOG(DEBUG) << "beginChild, now \"" << currentSolver_->name << "\" (" << currentSolver_
    << ") has " << currentSolver_->children.size() << "children: 0=" << currentSolver_->children[0];

  currentSolver_ = newChild;
}

//! indicate the end of the current child
void SolverStructureVisualizer::endChild()
{
  LOG(DEBUG) << "SolverStructureVisualizer::endChild() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  if (!currentSolver_)
    LOG(FATAL) << "endChild on invalid currentSolver";
  if (!currentSolver_->parent)
    LOG(FATAL) << "endChild, no parent";

  LOG(DEBUG) << "endChild, from \"" << currentSolver_->name << "\" back to \"" << currentSolver_->parent->name << "\".";


  currentSolver_ = currentSolver_->parent;
}

//! add the output connection information between two children to the current solver
void SolverStructureVisualizer::addOutputConnection(OutputConnection &outputConnection)
{
  LOG(DEBUG) << "SolverStructureVisualizer::addOutputConnection() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  // get output connector data from outputConnection
  const std::vector<OutputConnection::Connector> &connectorTerm1To2 = outputConnection.connectorTerm1To2();
  const std::vector<OutputConnection::Connector> &connectorTerm2To1 = outputConnection.connectorTerm2To1();

  // loop over connectors from term1 to term2
  for (int i = 0; i < connectorTerm1To2.size(); i++)
  {
    // if connector is not open
    if (connectorTerm1To2[i].index != -1)
    {
      solver_t::OutputConnection outputConnection;
      outputConnection.type = solver_t::OutputConnection::ab;
      outputConnection.fromSlot = i;
      outputConnection.toSlot = connectorTerm1To2[i].index;

      // if there exists a connection back between the same slots
      if (connectorTerm2To1.size() >= connectorTerm1To2[i].index)
      {
        if (connectorTerm2To1[connectorTerm1To2[i].index].index == i)
        {
          // if the connection type is to avoid copies in both directions
          if (connectorTerm1To2[i].avoidCopyIfPossible && connectorTerm2To1[connectorTerm1To2[i].index].avoidCopyIfPossible)
          {
            outputConnection.type = solver_t::OutputConnection::bidirectionalReuse;
          }
          else
          {
            outputConnection.type = solver_t::OutputConnection::bidirectionalCopy;
          }
        }
      }

      // add parsed output connection
      currentSolver_->outputConnections.push_back(outputConnection);
    }
  }

  // loop over connectors from term2 to term1
  for (int i = 0; i < connectorTerm2To1.size(); i++)
  {
    // if connector is not open
    if (connectorTerm2To1[i].index != -1)
    {
      solver_t::OutputConnection outputConnection;
      outputConnection.type = solver_t::OutputConnection::ba;
      outputConnection.fromSlot = i;
      outputConnection.toSlot = connectorTerm2To1[i].index;

      // if there exists a connection back between the same slots
      if (connectorTerm1To2.size() >= connectorTerm2To1[i].index)
      {
        if (connectorTerm1To2[connectorTerm2To1[i].index].index == i)
        {
          // if the connection type is to avoid copies in both directions
          if (connectorTerm2To1[i].avoidCopyIfPossible && connectorTerm1To2[connectorTerm2To1[i].index].avoidCopyIfPossible)
          {
            outputConnection.type = solver_t::OutputConnection::bidirectionalReuse;
          }
          else
          {
            outputConnection.type = solver_t::OutputConnection::bidirectionalCopy;
          }
        }
      }

      // add parsed output connection
      if (outputConnection.type != solver_t::OutputConnection::bidirectionalReuse
          && outputConnection.type != solver_t::OutputConnection::bidirectionalCopy)
      {
        currentSolver_->outputConnections.push_back(outputConnection);
      }
    }
  }
}

//! beginChild and addSolver will have no effect
void SolverStructureVisualizer::enable()
{
  nDisableCalls_--;
  if (nDisableCalls_ == 0)
  {
    enabled_ = true;
  }
  LOG(DEBUG) << "SolverStructureVisualizer::enable() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";
}

//! beginChild and addSolver will have no effect
void SolverStructureVisualizer::disable()
{
  enabled_ = false;
  nDisableCalls_++;

  LOG(DEBUG) << "SolverStructureVisualizer::disable() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";
}
