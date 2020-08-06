#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "slot_connection/slot_connector_data.h"

//! add a solver to the diagram
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void SolverStructureVisualizer::
setSlotConnectorData(std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>> slotConnectorData, bool isFromTuple)
{
  LOG(DEBUG) << "SolverStructureVisualizer::setSlotConnectorData() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";
  VLOG(1) << "at \"" << currentSolver_->name << "\" [" << currentSolver_ << "] setSlotConnectorData, isFromTuple=" << isFromTuple << ", n slots: " <<  currentSolver_->outputSlots.size();

  if (!enabled_)
    return;

  if (!slotConnectorData)
  {
    LOG(WARNING) << "In SolverStructureVisualizer::setSlotConnectorData, slotConnectorData is not set";
    return;
  }

  if (!isFromTuple)
  {
    LOG(DEBUG) << "at \"" << currentSolver_->name << "\" setSlotConnectorData, clear " << currentSolver_->outputSlots.size() << " slots";
    currentSolver_->outputSlots.clear();
  }

  // loop over ComponentOfFieldVariable entries as variable1 in slotConnectorData
  for (int i = 0; i < slotConnectorData->variable1.size(); i++)
  {
    Data::ComponentOfFieldVariable<FunctionSpaceType,nComponents1> &entry = slotConnectorData->variable1[i];

    // add the given field variable components to outputSlots
    solver_t::OutputSlot newOutputSlot;
    newOutputSlot.fieldVariableName = entry.values->name();
    newOutputSlot.componentName = entry.values->componentName(entry.componentNo);
    newOutputSlot.nComponents = nComponents1;
    newOutputSlot.variableNo = 1;
    newOutputSlot.meshDescription = entry.values->functionSpace()->getDescription();
    if (slotConnectorData->slotNames.size() > i)
      newOutputSlot.slotName = slotConnectorData->slotNames[i];

    currentSolver_->outputSlots.push_back(newOutputSlot);

    LOG(DEBUG) << "  slot " << entry.values->name() << "." << entry.values->componentName(entry.componentNo) << " (variable1), now " << currentSolver_->outputSlots.size();
  }

  // if the geometry is set, also add it to the list
  // do not do this, it is confusing
#if 0
  if (slotConnectorData->geometryField && !isFromTuple)
  {
    solver_t::OutputSlot newOutputSlot;
    newOutputSlot.fieldVariableName = "(geometry)";
    newOutputSlot.componentName = "(geometry)";
    newOutputSlot.nComponents = 3;
    newOutputSlot.variableNo = 1;
    newOutputSlot.meshDescription = slotConnectorData->geometryField->functionSpace()->getDescription();

    currentSolver_->outputSlots.push_back(newOutputSlot);
    LOG(DEBUG) << "  slot (geometry), now " << currentSolver_->outputSlots.size();
  }
#endif

  // loop over ComponentOfFieldVariable entries as variable2 in slotConnectorData
  for (int i = 0; i < slotConnectorData->variable2.size(); i++)
  {
    Data::ComponentOfFieldVariable<FunctionSpaceType,nComponents2> &entry = slotConnectorData->variable2[i];

    // add the given field variable components to outputSlots
    solver_t::OutputSlot newOutputSlot;
    newOutputSlot.fieldVariableName = entry.values->name();
    newOutputSlot.componentName = entry.values->componentName(entry.componentNo);
    newOutputSlot.nComponents = nComponents2;
    newOutputSlot.variableNo = 2;
    newOutputSlot.meshDescription = entry.values->functionSpace()->getDescription();

    if (slotConnectorData->slotNames.size() > slotConnectorData->variable1.size()+i)
      newOutputSlot.slotName = slotConnectorData->slotNames[slotConnectorData->variable1.size()+i];

    currentSolver_->outputSlots.push_back(newOutputSlot);

    LOG(DEBUG) << "  slot " << entry.values->name() << "." << entry.values->componentName(entry.componentNo) << " (variable2), now " << currentSolver_->outputSlots.size();
  }

  LOG(DEBUG) << "now there are " << currentSolver_->outputSlots.size() << " output slots";
}

template<typename T>
void SolverStructureVisualizer::
setSlotConnectorData(std::shared_ptr<std::vector<T>> slotConnectorData, bool isFromTuple)
{
  VLOG(1) << "at \"" << currentSolver_->name << "\" [" << currentSolver_ << "] setSlotConnectorData vector, isFromTuple=" << isFromTuple;
  if (!slotConnectorData)
  {
    LOG(WARNING) << "In SolverStructureVisualizer::setSlotConnectorData, slotConnectorData is not set";
    return;
  }

  // if slotConnectorData is a vector, only use the first entry
  if (!slotConnectorData->empty())
  {
    setSlotConnectorData((*slotConnectorData)[0], isFromTuple);
  }
}

template<typename SlotConnectorData1, typename SlotConnectorData2>
void SolverStructureVisualizer::
setSlotConnectorData(std::shared_ptr<std::tuple<SlotConnectorData1,SlotConnectorData2>> slotConnectorData, bool isFromTuple)
{
  VLOG(1) << "at \"" << currentSolver_->name << "\" [" << currentSolver_ << "] setSlotConnectorData tuple, isFromTuple=" << isFromTuple;
  setSlotConnectorData(std::get<0>(*slotConnectorData), true);
  setSlotConnectorData(std::get<1>(*slotConnectorData), true);
}
