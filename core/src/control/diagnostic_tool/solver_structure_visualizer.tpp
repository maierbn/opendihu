#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "data_management/output_connector_data.h"

//! add a solver to the diagram
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void SolverStructureVisualizer::
setOutputConnectorData(std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>> outputConnectorData, bool isFromTuple)
{
  LOG(DEBUG) << "SolverStructureVisualizer::setOutputConnectorData() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";
  VLOG(1) << "at \"" << currentSolver_->name << "\" [" << currentSolver_ << "] setOutputConnectorData, isFromTuple=" << isFromTuple << ", n slots: " <<  currentSolver_->outputSlots.size();

  if (!enabled_)
    return;

  if (!outputConnectorData)
  {
    LOG(WARNING) << "In SolverStructureVisualizer::setOutputConnectorData, outputConnectorData is not set";
    return;
  }

  if (!isFromTuple)
  {
    LOG(DEBUG) << "at \"" << currentSolver_->name << "\" setOutputConnectorData, clear " << currentSolver_->outputSlots.size() << " slots";
    currentSolver_->outputSlots.clear();
  }

  // loop over ComponentOfFieldVariable entries as variable1 in outputConnectorData
  for (int i = 0; i < outputConnectorData->variable1.size(); i++)
  {
    Data::ComponentOfFieldVariable<FunctionSpaceType,nComponents1> &entry = outputConnectorData->variable1[i];

    // add the given field variable components to outputSlots
    solver_t::OutputSlot newOutputSlot;
    newOutputSlot.fieldVariableName = entry.values->name();
    newOutputSlot.componentName = entry.values->componentName(entry.componentNo);
    newOutputSlot.nComponents = nComponents1;
    newOutputSlot.variableNo = 1;
    newOutputSlot.meshDescription = entry.values->functionSpace()->getDescription();

    currentSolver_->outputSlots.push_back(newOutputSlot);

    LOG(DEBUG) << "  slot " << entry.values->name() << "." << entry.values->componentName(entry.componentNo) << " (variable1), now " << currentSolver_->outputSlots.size();
  }

  // if the geometry is set, also add it to the list
  if (outputConnectorData->geometryField && !isFromTuple)
  {
    solver_t::OutputSlot newOutputSlot;
    newOutputSlot.fieldVariableName = "(geometry)";
    newOutputSlot.componentName = "(geometry)";
    newOutputSlot.nComponents = 3;
    newOutputSlot.variableNo = 1;
    newOutputSlot.meshDescription = outputConnectorData->geometryField->functionSpace()->getDescription();

    currentSolver_->outputSlots.push_back(newOutputSlot);
    LOG(DEBUG) << "  slot (geometry), now " << currentSolver_->outputSlots.size();
  }

  // loop over ComponentOfFieldVariable entries as variable2 in outputConnectorData
  for (int i = 0; i < outputConnectorData->variable2.size(); i++)
  {
    Data::ComponentOfFieldVariable<FunctionSpaceType,nComponents2> &entry = outputConnectorData->variable2[i];

    // add the given field variable components to outputSlots
    solver_t::OutputSlot newOutputSlot;
    newOutputSlot.fieldVariableName = entry.values->name();
    newOutputSlot.componentName = entry.values->componentName(entry.componentNo);
    newOutputSlot.nComponents = nComponents2;
    newOutputSlot.variableNo = 2;
    newOutputSlot.meshDescription = entry.values->functionSpace()->getDescription();

    currentSolver_->outputSlots.push_back(newOutputSlot);

    LOG(DEBUG) << "  slot " << entry.values->name() << "." << entry.values->componentName(entry.componentNo) << " (variable2), now " << currentSolver_->outputSlots.size();
  }

  LOG(DEBUG) << "now there are " << currentSolver_->outputSlots.size() << " output slots";
}

template<typename T>
void SolverStructureVisualizer::
setOutputConnectorData(std::shared_ptr<std::vector<T>> outputConnectorData, bool isFromTuple)
{
  VLOG(1) << "at \"" << currentSolver_->name << "\" [" << currentSolver_ << "] setOutputConnectorData vector, isFromTuple=" << isFromTuple;
  if (!outputConnectorData)
  {
    LOG(WARNING) << "In SolverStructureVisualizer::setOutputConnectorData, outputConnectorData is not set";
    return;
  }

  // if outputConnectorData is a vector, only use the first entry
  if (!outputConnectorData->empty())
  {
    setOutputConnectorData((*outputConnectorData)[0], isFromTuple);
  }
}

template<typename OutputConnectorData1, typename OutputConnectorData2>
void SolverStructureVisualizer::
setOutputConnectorData(std::shared_ptr<std::tuple<OutputConnectorData1,OutputConnectorData2>> outputConnectorData, bool isFromTuple)
{
  VLOG(1) << "at \"" << currentSolver_->name << "\" [" << currentSolver_ << "] setOutputConnectorData tuple, isFromTuple=" << isFromTuple;
  setOutputConnectorData(std::get<0>(*outputConnectorData), true);
  setOutputConnectorData(std::get<1>(*outputConnectorData), true);
}
