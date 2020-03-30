#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "data_management/output_connector_data.h"

//! add a solver to the diagram
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void SolverStructureVisualizer::
setOutputConnectorData(std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>> outputConnectorData)
{
  LOG(DEBUG) << "SolverStructureVisualizer::setOutputConnectorData() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  if (!outputConnectorData)
  {
    LOG(WARNING) << "In SolverStructureVisualizer::setOutputConnectorData, outputConnectorData is not set";
    return;
  }

  currentSolver_->outputSlots.clear();

  // loop over ComponentOfFieldVariable entries as variable1 in outputConnectorData
  for (int i = 0; i < outputConnectorData->variable1.size(); i++)
  {
    Data::ComponentOfFieldVariable<FunctionSpaceType,nComponents1> &entry = outputConnectorData->variable1[i];

    // add the given field variable components to outputSlots
    currentSolver_->outputSlots.emplace_back();
    currentSolver_->outputSlots.back().fieldVariableName = entry.values->name();
    currentSolver_->outputSlots.back().componentName = entry.values->componentName(entry.componentNo);
    currentSolver_->outputSlots.back().nComponents = nComponents1;
    currentSolver_->outputSlots.back().variableNo = 1;
  }

  // loop over ComponentOfFieldVariable entries as variable2 in outputConnectorData
  for (int i = 0; i < outputConnectorData->variable2.size(); i++)
  {
    Data::ComponentOfFieldVariable<FunctionSpaceType,nComponents2> &entry = outputConnectorData->variable2[i];

    // add the given field variable components to outputSlots
    currentSolver_->outputSlots.emplace_back();
    currentSolver_->outputSlots.back().fieldVariableName = entry.values->name();
    currentSolver_->outputSlots.back().componentName = entry.values->componentName(entry.componentNo);
    currentSolver_->outputSlots.back().nComponents = nComponents2;
    currentSolver_->outputSlots.back().variableNo = 2;
  }

  LOG(DEBUG) << "added " << currentSolver_->outputSlots.size() << " output slots";
}

template<typename T>
void SolverStructureVisualizer::
setOutputConnectorData(std::shared_ptr<std::vector<T>> outputConnectorData)
{
  if (!outputConnectorData)
  {
    LOG(WARNING) << "In SolverStructureVisualizer::setOutputConnectorData, outputConnectorData is not set";
    return;
  }

  // if outputConnectorData is a vector, only use the first entry
  if (!outputConnectorData->empty())
  {
    setOutputConnectorData((*outputConnectorData)[0]);
  }
}
