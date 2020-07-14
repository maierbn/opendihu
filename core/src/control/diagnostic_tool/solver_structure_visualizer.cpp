#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "output_connector_data_transfer/output_connector_data.h"
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

void SolverStructureVisualizer::addSolver(std::string name, bool hasInternalConnectionToFirstNestedSolver, bool hasInternalConnectionToSecondNestedSolver)
{
  LOG(DEBUG) << "SolverStructureVisualizer::addSolver(\"" <<  name << "\",hasInternalConnectionToFirstNestedSolver=" << hasInternalConnectionToFirstNestedSolver 
    << "," << hasInternalConnectionToSecondNestedSolver << ") under \"" << currentSolver_->parent->name << "\", nDisableCalls_: " << nDisableCalls_
    << ", enabled: " << enabled_  << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  if (!currentSolver_)
    LOG(FATAL) << "addSolver on invalid currentSolver";

  if (currentSolver_->name != "" && currentSolver_->name != name)
  {
    LOG(WARNING) << "SolverStructureVisualizer::addSolver(\"" <<  name << "\") under \"" << currentSolver_->parent->name << "\", nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_
      << " overwrites solver name \"" << currentSolver_->name << "\".";
  }
  currentSolver_->name = name;
  currentSolver_->hasInternalConnectionToFirstNestedSolver = hasInternalConnectionToFirstNestedSolver;
  currentSolver_->hasInternalConnectionToSecondNestedSolver = hasInternalConnectionToSecondNestedSolver;

  LOG(DEBUG) << "addSolver \"" << name << "\" [" << currentSolver_ << "].";
}

//! add a description for the current that will be included in the visualization, to be called after addSolver
void SolverStructureVisualizer::setSolverDescription(std::string description)
{
  LOG(DEBUG) << "SolverStructureVisualizer::setSolverDescription(\"" <<  description << "\")"
    << " under \"" << currentSolver_->parent->name << "\", nDisableCalls_: " << nDisableCalls_
    << ", enabled: " << enabled_  << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  if (!currentSolver_)
    LOG(FATAL) << "addSolver on invalid currentSolver";

  currentSolver_->description = description;
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

  LOG(DEBUG) << "beginChild, now \"" << currentSolver_->name << "\" [" << currentSolver_
    << "] has " << currentSolver_->children.size() << "children: 0=[" << currentSolver_->children[0] << "]";

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

  LOG(DEBUG) << "endChild, from \"" << currentSolver_->name << "\"  [" << currentSolver_ << "] back to \"" << currentSolver_->parent->name << "\" [" << currentSolver_->parent << "].";


  currentSolver_ = currentSolver_->parent;
}

//! add the output connection information between two children to the current solver
void SolverStructureVisualizer::addOutputConnection(std::shared_ptr<OutputConnection> outputConnection)
{
  if (currentSolver_)
    LOG(DEBUG) << "SolverStructureVisualizer::addOutputConnection() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  currentSolver_->outputConnection = outputConnection;
}

void SolverStructureVisualizer::parseOutputConnection(std::shared_ptr<solver_t> currentSolver)
{
  VLOG(1) << "parseOutputConnection";

  if (!currentSolver->outputConnection)
  {
    VLOG(1) << "output connection not set";
    return;
  }

  LOG(DEBUG) << "currentSolver->outputConnection: " << currentSolver->outputConnection;

  currentSolver->outputConnections.clear();

  // get output connector data from outputConnection
  const std::vector<OutputConnection::Connector> &connectorTerm1To2 = currentSolver->outputConnection->connectorForVisualizerTerm1To2();
  const std::vector<OutputConnection::Connector> &connectorTerm2To1 = currentSolver->outputConnection->connectorForVisualizerTerm2To1();

  // loop over connectors from term1 to term2
  for (int i = 0; i < connectorTerm1To2.size(); i++)
  {
    LOG(DEBUG) << "term1 -> term2 i=" << i << ", map to " << connectorTerm1To2[i].index << ", avoidCopyIfPossible: " << connectorTerm1To2[i].avoidCopyIfPossible;

    // if connector is not open
    if (connectorTerm1To2[i].index != -1)
    {
      solver_t::OutputConnectionRepresentation outputConnection;
      outputConnection.type = solver_t::OutputConnectionRepresentation::ab;
      outputConnection.fromSlot = i;
      outputConnection.toSlot = connectorTerm1To2[i].index;

      // if there exists a connection back between the same slots
      if (connectorTerm2To1.size() > outputConnection.toSlot)
      {
        if (connectorTerm2To1[outputConnection.toSlot].index == i)
        {
          // if the connection type is to avoid copies in both directions
          if (connectorTerm1To2[i].avoidCopyIfPossible && connectorTerm2To1[outputConnection.toSlot].avoidCopyIfPossible)
          {
            outputConnection.type = solver_t::OutputConnectionRepresentation::bidirectionalReuse;
          }
          else
          {
            outputConnection.type = solver_t::OutputConnectionRepresentation::bidirectionalCopy;
          }
        }
      }

      // add parsed output connection
      currentSolver->outputConnections.push_back(outputConnection);
    }
  }

  // loop over connectors from term2 to term1
  for (int i = 0; i < connectorTerm2To1.size(); i++)
  {
    LOG(DEBUG) << "term2 -> term1 i=" << i << ", map to " << connectorTerm2To1[i].index << ", avoidCopyIfPossible: " << connectorTerm2To1[i].avoidCopyIfPossible;

    // if connector is not open
    if (connectorTerm2To1[i].index != -1)
    {
      solver_t::OutputConnectionRepresentation outputConnection;
      outputConnection.type = solver_t::OutputConnectionRepresentation::ba;
      outputConnection.fromSlot = i;
      outputConnection.toSlot = connectorTerm2To1[i].index;

      // if there exists a connection back between the same slots
      if (connectorTerm1To2.size() > outputConnection.toSlot)
      {
        if (connectorTerm1To2[outputConnection.toSlot].index == i)
        {
          // if the connection type is to avoid copies in both directions
          if (connectorTerm2To1[i].avoidCopyIfPossible && connectorTerm1To2[outputConnection.toSlot].avoidCopyIfPossible)
          {
            outputConnection.type = solver_t::OutputConnectionRepresentation::bidirectionalReuse;
          }
          else
          {
            outputConnection.type = solver_t::OutputConnectionRepresentation::bidirectionalCopy;
          }
        }
      }

      // add parsed output connection
      if (outputConnection.type != solver_t::OutputConnectionRepresentation::bidirectionalReuse
          && outputConnection.type != solver_t::OutputConnectionRepresentation::bidirectionalCopy)
      {
        currentSolver->outputConnections.push_back(outputConnection);
      }
    }
  }
}

//! add connections between slots that occur within the same solver
void SolverStructureVisualizer::addSlotMapping(int slotNoFrom, int slotNoTo)
{
  if (!enabled_)
    return;

  // check if this mapping does not yet exist
  bool mappingExists = false;
  for (const solver_t::OutputConnectionRepresentation &mappingWithinSolver : currentSolver_->mappingsWithinSolver)
  {
    if (mappingWithinSolver.fromSlot == slotNoFrom && mappingWithinSolver.toSlot == slotNoTo)
    {
      mappingExists = true;
      break;
    }
  }

  if (mappingExists)
    return;

  solver_t::OutputConnectionRepresentation mappingWithinSolver;
  mappingWithinSolver.fromSlot = slotNoFrom;
  mappingWithinSolver.toSlot = slotNoTo;
  currentSolver_->mappingsWithinSolver.push_back(mappingWithinSolver);
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

//! produce the resulting file
std::string SolverStructureVisualizer::
getDiagram()
{
  if (!enabled_)
    return std::string("");

  std::string diagram;

  // only produce file on rank 0
  if (DihuContext::ownRankNoCommWorld() == 0)
  {

    // collect all information

    // set to root of nested solvers to start printing from there
    DiagramGenerator diagramGenerator;
    diagramGenerator.initialize(solverRoot_);

    if (solverRoot_)
      LOG(DEBUG) << "SolverStructureVisualizer::getDiagram() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << solverRoot_->name << "\".";

    // call generateDiagram on the nested solvers
    diagram = diagramGenerator.generateFinalDiagram();

  }

  return diagram;
}
