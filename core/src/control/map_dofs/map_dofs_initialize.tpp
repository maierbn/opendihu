#include "control/map_dofs/map_dofs.h"

#include "slot_connection/data_helper/slot_connector_data_helper.h"

namespace Control
{

//! constructor, gets the DihuContext object which contains all python settings
template<typename FunctionSpaceType, typename NestedSolverType>
MapDofs<FunctionSpaceType,NestedSolverType>::
MapDofs(DihuContext context) :
  context_(context["MapDofs"]), nestedSolver_(context_), data_(context_)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();
}

//! initialize field variables and everything needed for the dofs mapping
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
initialize()
{
  LOG(DEBUG) << "initialize MapDofs";

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("MapDofs", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // parse description for solverStructureVisualizer, if there was any
  std::string description;
  if (this->specificSettings_.hasKey("description"))
    description = this->specificSettings_.getOptionString("description", "");
  DihuContext::solverStructureVisualizer()->setSolverDescription(description);

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // initialize solver
  nestedSolver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // parse settings
  int nAdditionalFieldVariables = this->specificSettings_.getOptionInt("nAdditionalFieldVariables", 0, PythonUtility::NonNegative);

  // create function space to use for the additional field variables
  std::shared_ptr<FunctionSpaceType> functionSpace = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);
  data_.setFunctionSpace(functionSpace);

  // initialize data, i.e. create additional field variables
  data_.initialize(nAdditionalFieldVariables, nestedSolver_);

  // get all available slot names
  slotNames_.clear();
  SlotConnectorDataHelper<SlotConnectorDataType>::getSlotNames(getSlotConnectorData(), slotNames_);

  parseMappingFromSettings("beforeComputation", mappingsBeforeComputation_);
  parseMappingFromSettings("afterComputation", mappingsAfterComputation_);

  // prepare communication by initializing mappings
  initializeCommunication(mappingsBeforeComputation_);
  initializeCommunication(mappingsAfterComputation_);

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());

  // add mappings to be visualized by solverStructureVisualizer
  using connection_t = SolverStructureVisualizer::solver_t::SlotsConnectionRepresentation;
  for (const DofsMappingType &mapping : mappingsBeforeComputation_)
  {
    for (int connectorSlotNoTo : mapping.connectorSlotNosTo)
    {
      DihuContext::solverStructureVisualizer()->addSlotMapping(mapping.connectorSlotNoFrom, connectorSlotNoTo, connection_t::internalBeforeComputation);
    }
  }
  for (const DofsMappingType &mapping : mappingsAfterComputation_)
  {
    for (int connectorSlotNoTo : mapping.connectorSlotNosTo)
    {
      DihuContext::solverStructureVisualizer()->addSlotMapping(mapping.connectorSlotNoFrom, connectorSlotNoTo, connection_t::internalAfterComputation);
    }
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
parseMappingFromSettings(std::string settingsKey, std::vector<DofsMappingType> &mappings)
{

  // parse list items under settingsKey="beforeComputation" or "afterComputation"
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  if (listPy == Py_None)
    return;

  std::vector<PyObject *> list = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);

  PythonConfig mappingSpecification(this->specificSettings_, settingsKey);

  // loop over items of the list under "beforeComputation" or "afterComputation"
  for (int i = 0; i < list.size(); i++)
  {
    PythonConfig currentMappingSpecification(mappingSpecification, i);

    DofsMappingType newDofsMapping;

    // output errors if the old options are used
    if (currentMappingSpecification.hasKey("fromOutputConnectorSlotNo"))
    {
      LOG(ERROR) << currentMappingSpecification << ": Option \"fromOutputConnectorSlotNo\" has been renamed to \"fromConnectorSlot\".";
    }
    if (currentMappingSpecification.hasKey("toOutputConnectorSlotNo"))
    {
      LOG(ERROR) << currentMappingSpecification << ": Option \"toOutputConnectorSlotNo\" has been renamed to \"toConnectorSlots\".";
    }
    if (currentMappingSpecification.hasKey("fromOutputConnectorArrayIndex"))
    {
      LOG(ERROR) << currentMappingSpecification << ": Option \"fromOutputConnectorArrayIndex\" has been renamed to \"fromSlotConnectorArrayIndex\".";
    }
    if (currentMappingSpecification.hasKey("toOutputConnectorArrayIndex"))
    {
      LOG(ERROR) << currentMappingSpecification << ": Option \"toOutputConnectorArrayIndex\" has been renamed to \"toSlotConnectorArrayIndex\".";
    }

    // parse options of this item
    // parse fromConnectorSlot, either as integer value (then it is the no.) or as string (then it is the slot name)
    PyObject *fromConnectorSlotPy = currentMappingSpecification.getOptionPyObject("fromConnectorSlot");
    if (PyLong_Check(fromConnectorSlotPy))
    {
      // parse slot no
      newDofsMapping.connectorSlotNoFrom = currentMappingSpecification.getOptionInt("fromConnectorSlot", 0, PythonUtility::NonNegative);
    }
    else
    {
      // if the option is not an integer, it is interpreted as string as the slot name
      std::string slotName = currentMappingSpecification.getOptionString("fromConnectorSlot", "");

      // check if this slot name exists
      std::vector<std::string>::iterator iter = std::find(slotNames_.begin(), slotNames_.end(), slotName);
      if (iter == slotNames_.end())
      {
        LOG(FATAL) << currentMappingSpecification << "[\"fromConnectorSlot\"] = " << slotName << " is interpreted as slot name, but there is no such slot. Available slot names: " << slotNames_;
      }
      newDofsMapping.connectorSlotNoFrom = std::distance(slotNames_.begin(), iter);
    }

    // parse toConnectorSlots, either as integer value (then it is the no.) or as string (then it is the slot name)
    PyObject *toConnectorSlotsPy = currentMappingSpecification.getOptionPyObject("toConnectorSlots");
    if (PyList_Check(toConnectorSlotsPy))
    {
      int nEntries = PyList_Size(toConnectorSlotsPy);
      if (nEntries == 0)
        LOG(FATAL) << currentMappingSpecification << "[\"toConnectorSlots\"] is an empty list.";

      std::vector<PyObject *> toConnectorSlotsItemsPy;
      currentMappingSpecification.getOptionVector<PyObject *>("toConnectorSlots", toConnectorSlotsItemsPy);

      for (PyObject *toConnectorSlotPy : toConnectorSlotsItemsPy)
      {
        if (PyLong_Check(toConnectorSlotPy))
        {
          // parse slot no
          newDofsMapping.connectorSlotNosTo.push_back(PythonUtility::convertFromPython<int>::get(toConnectorSlotPy, 0));
        }
        else
        {
          // if the option is not an integer, it is interpreted as string as the slot name
          std::string slotName = PythonUtility::convertFromPython<std::string>::get(toConnectorSlotPy);

          // check if this slot name exists
          std::vector<std::string>::iterator iter = std::find(slotNames_.begin(), slotNames_.end(), slotName);
          if (iter == slotNames_.end())
          {
            LOG(FATAL) << currentMappingSpecification << "[\"toConnectorSlots\"][...] = " << slotName << " is interpreted as slot name, but there is no such slot. Available slot names: " << slotNames_;
          }
          newDofsMapping.connectorSlotNosTo.push_back(std::distance(slotNames_.begin(), iter));
        }
      }
    }
    else if (PyLong_Check(toConnectorSlotsPy))
    {
      // parse slot no
      newDofsMapping.connectorSlotNosTo.push_back(currentMappingSpecification.getOptionInt("toConnectorSlots", 0, PythonUtility::NonNegative));
    }
    else
    {
      // if the option is not an integer, it is interpreted as string as the slot name
      std::string slotName = currentMappingSpecification.getOptionString("toConnectorSlots", "");

      // check if this slot name exists
      std::vector<std::string>::iterator iter = std::find(slotNames_.begin(), slotNames_.end(), slotName);
      if (iter == slotNames_.end())
      {
        LOG(FATAL) << currentMappingSpecification << "[\"toConnectorSlots\"] = " << slotName << " is interpreted as slot name, but there is no such slot. Available slot names: " << slotNames_;
      }
      newDofsMapping.connectorSlotNosTo.push_back(std::distance(slotNames_.begin(), iter));
    }

    newDofsMapping.slotConnectorArrayIndexFrom = currentMappingSpecification.getOptionInt("fromSlotConnectorArrayIndex", 0, PythonUtility::NonNegative);
    newDofsMapping.slotConnectorArrayIndexTo   = currentMappingSpecification.getOptionInt("toSlotConnectorArrayIndex",   0, PythonUtility::NonNegative);

    // parse if dof nos are given as global or local nos
    std::string fromDofNosNumbering = currentMappingSpecification.getOptionString("fromDofNosNumbering", "global");
    std::string toDofNosNumbering   = currentMappingSpecification.getOptionString("toDofNosNumbering",   "local");

    // parse dofNoIsGlobalFrom
    if (fromDofNosNumbering == "global")
    {
      newDofsMapping.dofNoIsGlobalFrom = true;
    }
    else if (fromDofNosNumbering == "local")
    {
      newDofsMapping.dofNoIsGlobalFrom = false;
    }
    else
    {
      LOG(ERROR) << currentMappingSpecification << "[\"fromDofNosNumbering\"] is \"" << fromDofNosNumbering << ", but has to be \"global\" or \"local\". Assuming \"local\".";
      newDofsMapping.dofNoIsGlobalFrom = false;
    }

    // parse dofNoIsGlobalTo
    if (toDofNosNumbering == "global")
    {
      newDofsMapping.dofNoIsGlobalTo = true;
    }
    else if (toDofNosNumbering == "local")
    {
      newDofsMapping.dofNoIsGlobalTo = false;
    }
    else
    {
      LOG(ERROR) << currentMappingSpecification << "[\"toDofNosNumbering\"] is \"" << toDofNosNumbering << ", but has to be \"global\" or \"local\". Assuming \"local\".";
      newDofsMapping.dofNoIsGlobalTo = false;
    }

    // parse mode
    std::string mode = currentMappingSpecification.getOptionString("mode", "copyLocal");
    if (mode == "copyLocal")
    {
      newDofsMapping.mode = DofsMappingType::modeCopyLocal;
    }
    else if (mode == "copyLocalIfPositive")
    {
      newDofsMapping.mode = DofsMappingType::modeCopyLocalIfPositive;
    }
    else if (mode == "communicate")
    {
      newDofsMapping.mode = DofsMappingType::modeCommunicate;
    }
    else if (mode == "localSetIfAboveThreshold")
    {
      newDofsMapping.mode = DofsMappingType::modeLocalSetIfAboveThreshold;
      newDofsMapping.thresholdValue = currentMappingSpecification.getOptionDouble("thresholdValue", 0);
      newDofsMapping.valueToSet = currentMappingSpecification.getOptionDouble("valueToSet", 0);
    }
    else if (mode == "callback")
    {
      newDofsMapping.mode = DofsMappingType::modeCallback;

      // parse callback function
      newDofsMapping.callback = currentMappingSpecification.getOptionPyObject("callback");

      // parse input and output dof nos
      currentMappingSpecification.getOptionVector<dof_no_t>("inputDofs", newDofsMapping.inputDofs);
      PyObject *outputDofsObject = currentMappingSpecification.getOptionPyObject("outputDofs");
      newDofsMapping.outputDofs = PythonUtility::convertFromPython<std::vector<std::vector<dof_no_t>>>::get(outputDofsObject);

      // check that the number of output dof lists matches the number of connector slots
      int nOutputDofSlots = newDofsMapping.outputDofs.size();
      int nConnectorSlotNos = newDofsMapping.connectorSlotNosTo.size();
      if (nOutputDofSlots != nConnectorSlotNos)
      {
        LOG(FATAL) << currentMappingSpecification << "[\"toConnectorSlots\"] specifies " << nConnectorSlotNos << " connector slots, but "
         << currentMappingSpecification << "[\"outputDofs\"] contains " << nOutputDofSlots << " lists of dofs.";
      }

      // create python list [fromSlotNo, [toSlotNos], fromArrayIndex, toArrayIndex]
      // this list is needed as argument for the callback function
      std::vector<int> slotNosList{
        newDofsMapping.connectorSlotNoFrom, newDofsMapping.connectorSlotNosTo[0],
        newDofsMapping.slotConnectorArrayIndexFrom, newDofsMapping.slotConnectorArrayIndexTo
      };
      newDofsMapping.slotNosPy = PythonUtility::convertToPython<std::vector<int>>::get(slotNosList);

      // set the second list item in [fromSlotNo, [...], fromArrayIndex, toArrayIndex] to [toSlotNo0, toSlotNo1, ...]
      std::vector<int> connectorSlotNosToList = newDofsMapping.connectorSlotNosTo;
      PyObject *connectorSlotNosToListPy = PythonUtility::convertToPython<std::vector<int>>::get(connectorSlotNosToList);

      PyList_SetItem(newDofsMapping.slotNosPy, Py_ssize_t(1), connectorSlotNosToListPy);

      // initialize the user buffer
      newDofsMapping.buffer = PyDict_New();

      // initialize outputValues as [[None, None, ..., None],[...]]
      newDofsMapping.outputValuesPy = PyList_New(Py_ssize_t(nOutputDofSlots));

      for (int outputDofSlotIndex = 0; outputDofSlotIndex < nOutputDofSlots; outputDofSlotIndex++)
      {
        int nInitialOutputValues = newDofsMapping.outputDofs[outputDofSlotIndex].size();
        PyObject *slotList = PyList_New(Py_ssize_t(nInitialOutputValues));

        for (int i = 0; i < nInitialOutputValues; i++)
          PyList_SetItem(slotList, Py_ssize_t(i), Py_None);

        PyList_SetItem(newDofsMapping.outputValuesPy, Py_ssize_t(outputDofSlotIndex), slotList);
      }
    }
    else
    {
      LOG(FATAL) << currentMappingSpecification << "[\"mode\"] is \"" << mode << "\", but allowed values are "
        <<"\"copyLocal\", \"copyLocalIfPositive\", \"callback\" and \"communicate\".";
    }

    PyObject *object = currentMappingSpecification.getOptionPyObject("dofsMapping");
    if (object != Py_None)
      newDofsMapping.dofsMapping = PythonUtility::convertFromPython<std::map<int,std::vector<int>>>::get(object);

    LOG(DEBUG) << "Parsed mapping slots " << newDofsMapping.connectorSlotNoFrom << " (arrayIndex " << newDofsMapping.slotConnectorArrayIndexFrom
      << ") -> " << newDofsMapping.connectorSlotNosTo << " (arrayIndex " << newDofsMapping.slotConnectorArrayIndexTo << "), dofNos "
      << (newDofsMapping.dofNoIsGlobalFrom? "global" : "local") << " -> " << (newDofsMapping.dofNoIsGlobalTo? "global" : "local") << ", mode: " << mode;
    LOG(DEBUG) << "dofsMapping: " << newDofsMapping.dofsMapping;

    mappings.push_back(newDofsMapping);
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
std::shared_ptr<Partition::MeshPartitionBase> MapDofs<FunctionSpaceType,NestedSolverType>::
getMeshPartitionBase(int slotNo, int arrayIndex)
{
  /*
  typedef std::tuple<
    std::shared_ptr<typename NestedSolverType::SlotConnectorDataType>,
    std::shared_ptr<SlotConnectorData<FunctionSpaceType,1>>
  > SlotConnectorDataType;
  */

  // need function space of affected field variables
  int nSlotsNestedSolver = SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::nSlots(
    std::get<0>(*data_.getSlotConnectorData())
  );

  if (slotNo < nSlotsNestedSolver)
  {
    return SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::getMeshPartitionBase(
      std::get<0>(*data_.getSlotConnectorData()), slotNo, arrayIndex
    );
  }
  else
  {
    int index = slotNo - nSlotsNestedSolver;
    int nFieldVariables = std::get<1>(*data_.getSlotConnectorData())->variable1.size();
    if (index < nFieldVariables)
    {

      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
        = std::get<1>(*data_.getSlotConnectorData())->variable1[index].values;

      if (fieldVariable->functionSpace())
      {
        return fieldVariable->functionSpace()->meshPartitionBase();
      }
      else
      {
        LOG(ERROR) << this->specificSettings_ << " Slot " << slotNo << " refers to a not initialized field variable.";
      }
    }
    else
    {
      LOG(ERROR) << this->specificSettings_ << " Slot " << slotNo << " refers to a non-existing field variable, "
        << "index=" << index << " but number of available field variables: " << nFieldVariables;
    }
  }

  return nullptr;
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
slotGetValues(int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values)
{
  //LOG(DEBUG) << "   slotGetValues(slotNo=" << slotNo << ", arrayIndex=" << arrayIndex << ", " << dofNosLocal.size() << " dofs: " << dofNosLocal;
  
  // need function space of affected field variables
  int nSlotsNestedSolver = SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::nSlots(
    std::get<0>(*data_.getSlotConnectorData())
  );
  if (slotNo < nSlotsNestedSolver)
  {
    SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::slotGetValues(
      std::get<0>(*data_.getSlotConnectorData()), slotNo, arrayIndex, dofNosLocal, values
    );
    //LOG(DEBUG) << "   got values from nested solver slots: " << values;
  }
  else
  {
    // get the field variable and component no for this slot
    int index = slotNo - nSlotsNestedSolver;

    int nAdditionalSlots = std::get<1>(*data_.getSlotConnectorData())->variable1.size();
    if (index < 0 || index >= nAdditionalSlots)
    {
      LOG(FATAL) << "In MapDofs, cannot get values for slot " << slotNo << ", number of slots: " << nSlotsNestedSolver + nAdditionalSlots
        << " (nested solver has " << nSlotsNestedSolver << " slots and there " << (nAdditionalSlots==1? "is ": "are ") << nAdditionalSlots << " additional slot" << (nAdditionalSlots==1? "" : "s") << ")";
    }

    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
      = std::get<1>(*data_.getSlotConnectorData())->variable1[index].values;

    int componentNo = std::get<1>(*data_.getSlotConnectorData())->variable1[index].componentNo;  // should be 0

    if (!fieldVariable)
    {
      LOG(FATAL) << "In MapDofs, cannot get values, field variable is not set for slot " << slotNo << ", number of slots: " << nSlotsNestedSolver + nAdditionalSlots
        << " (nested solver has " << nSlotsNestedSolver << " slots and there " << (nAdditionalSlots==1? "is ": "are ") << nAdditionalSlots << " additional slot" << (nAdditionalSlots==1? "" : "s") << ")";
    }

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    //LOG(DEBUG) << "   got values from additional field variable no " << index << ", \"" << fieldVariable->name() << "\", component " << componentNo
    //  << ", got values " << values;
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
bool MapDofs<FunctionSpaceType,NestedSolverType>::
slotSetValues(int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode)
{
  LOG(DEBUG) << "   slotSetValues(slotNo=" << slotNo << ", arrayIndex=" << arrayIndex << ", " << dofNosLocal.size() << " dofs: " << dofNosLocal;// << ", values: " << values;
  
  // need function space of affected field variables
  int nSlotsNestedSolver = SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::nSlots(
    std::get<0>(*data_.getSlotConnectorData())
  );
  if (slotNo < nSlotsNestedSolver)
  {
    return SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::slotSetValues(
      std::get<0>(*data_.getSlotConnectorData()), slotNo, arrayIndex, dofNosLocal, values, petscInsertMode
    );
    //LOG(DEBUG) << "   set values in nested solver slots";
  }
  else
  {
    // get the field variable and component no for this slot
    int index = slotNo - nSlotsNestedSolver;

    int nAdditionalSlots = std::get<1>(*data_.getSlotConnectorData())->variable1.size();
    if (index < 0 || index >= nAdditionalSlots)
    {
      LOG(FATAL) << "In MapDofs, cannot set values for slot " << slotNo << ", number of slots: " << nSlotsNestedSolver + nAdditionalSlots
        << " (nested solver has " << nSlotsNestedSolver << " slots and there " << (nAdditionalSlots==1? "is ": "are ") << nAdditionalSlots << " additional slot" << (nAdditionalSlots==1? "" : "s") << ")";
    }

    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
      = std::get<1>(*data_.getSlotConnectorData())->variable1[index].values;

    int componentNo = std::get<1>(*data_.getSlotConnectorData())->variable1[index].componentNo;  // should be 0

    if (!fieldVariable)
    {
      LOG(FATAL) << "In MapDofs, cannot set values, field variable is not set for slot " << slotNo << ", number of slots: " << nSlotsNestedSolver + nAdditionalSlots
        << " (nested solver has " << nSlotsNestedSolver << " slots and there " << (nAdditionalSlots==1? "is ": "are ") << nAdditionalSlots << " additional slot" << (nAdditionalSlots==1? "" : "s") << ")";
    }

    LOG(DEBUG) << "   set values in additional fieldVariable no " << index << ", \"" << fieldVariable->name() << "\", component " << componentNo
      << " on mesh \"" << fieldVariable->functionSpace()->meshName() << "\", set dofs " << dofNosLocal << " to values " << values;
      
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    return true;
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
slotSetRepresentationGlobal(int slotNo, int arrayIndex)
{
  // need function space of affected field variables
  int nSlotsNestedSolver = SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::nSlots(
    std::get<0>(*data_.getSlotConnectorData())
  );
  if (slotNo < nSlotsNestedSolver)
  {
    SlotConnectorDataHelper<typename NestedSolverType::SlotConnectorDataType>::slotSetRepresentationGlobal(
      std::get<0>(*data_.getSlotConnectorData()), slotNo, arrayIndex
    );
  }
  else
  {
    // get the field variable and component no for this slot
    int index = slotNo - nSlotsNestedSolver;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
      = std::get<1>(*data_.getSlotConnectorData())->variable1[index].values;
    fieldVariable->setRepresentationGlobal();
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
initializeCommunication(std::vector<DofsMappingType> &mappings)
{
  int ownRankNo = data_.functionSpace()->meshPartition()->rankSubset()->ownRankNo();

  // loop over mappings
  for (DofsMappingType &mapping : mappings)
  {
    // get the mesh partition of the "from" field variable
    std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBaseFrom = getMeshPartitionBase(mapping.connectorSlotNoFrom, mapping.slotConnectorArrayIndexFrom);

    if (!meshPartitionBaseFrom)
      LOG(FATAL) << "MapDofs: Could not get mesh partition for \"from\" function space for mapping " << mapping.connectorSlotNoFrom
        << " -> " << mapping.connectorSlotNosTo;

    if (mapping.mode == DofsMappingType::modeCallback)
    {
      if (mapping.dofNoIsGlobalFrom)
      {
        // transform the global input dofs to local input dofs
        std::vector<dof_no_t> localDofsMapping;
        for (std::vector<dof_no_t>::iterator iter = mapping.inputDofs.begin(); iter != mapping.inputDofs.end(); iter++)
        {
          global_no_t dofNoGlobalPetsc = *iter;
          bool isLocal;
          dof_no_t dofNoLocal = meshPartitionBaseFrom->getDofNoLocal(dofNoGlobalPetsc, isLocal);

          if (isLocal)
          {
            localDofsMapping.push_back(dofNoLocal);
          }
        }

        LOG(DEBUG) << "inputDofs global: " << mapping.inputDofs << ", transformed to local: " << localDofsMapping;

        // assign to dofs mapping
        mapping.inputDofs = localDofsMapping;
      }

      if (mapping.dofNoIsGlobalTo)
      {
        // transform the global input dofs to local input dofs
        std::vector<std::vector<dof_no_t>> localDofsMapping;
        localDofsMapping.resize(mapping.connectorSlotNosTo.size());

        // loop over the "to" slots
        for (int toSlotIndex = 0; toSlotIndex < mapping.connectorSlotNosTo.size(); toSlotIndex++)
        {
          int connectorSlotNoTo = mapping.connectorSlotNosTo[toSlotIndex];

          // get the mesh partition of the "to" field variable
          std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBaseTo = getMeshPartitionBase(connectorSlotNoTo, mapping.slotConnectorArrayIndexTo);

          if (!meshPartitionBaseTo)
            LOG(FATAL) << "MapDofs: Could not get mesh partition for \"to\" function space for mapping " << mapping.connectorSlotNoFrom
              << " -> " << connectorSlotNoTo;

          for (std::vector<dof_no_t>::iterator iter = mapping.outputDofs[toSlotIndex].begin(); iter != mapping.outputDofs[toSlotIndex].end(); iter++)
          {
            global_no_t dofNoGlobalPetsc = *iter;
            bool isLocal;
            dof_no_t dofNoLocal = meshPartitionBaseTo->getDofNoLocal(dofNoGlobalPetsc, isLocal);

            if (isLocal)
            {
              localDofsMapping[toSlotIndex].push_back(dofNoLocal);
            }
          }

          LOG(DEBUG) << "outputDofs global: " << mapping.outputDofs[toSlotIndex] << ", transformed to local: " << localDofsMapping;
        }
        // assign to dofs mapping
        mapping.outputDofs = localDofsMapping;
      }
      LOG(DEBUG) << "for modeCallback, initialized inputDofs: " << mapping.inputDofs << ", outputDofs: " << mapping.outputDofs;
    }
    else
    {
      // mode is not callback, this means we have exactly one output slot

      // get the mesh partition of the "to" field variable
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBaseTo = getMeshPartitionBase(mapping.connectorSlotNosTo[0], mapping.slotConnectorArrayIndexTo);

      if (!meshPartitionBaseTo)
        LOG(FATAL) << "MapDofs: Could not get mesh partition for \"to\" function space for mapping " << mapping.connectorSlotNoFrom
          << " -> " << mapping.connectorSlotNosTo[0];

      // convert from dof nos from global to local nos, discard values that are not on the local subdomain
      if (mapping.dofNoIsGlobalFrom)
      {
        // construct a new dofsMapping where all keys are local values and only include those that are local
        std::map<int,std::vector<int>> localDofsMapping;
        for (std::map<int,std::vector<int>>::iterator iter = mapping.dofsMapping.begin(); iter != mapping.dofsMapping.end(); iter++)
        {
          global_no_t dofNoGlobalPetsc = iter->first;
          bool isLocal;
          dof_no_t dofNoLocal = meshPartitionBaseFrom->getDofNoLocal(dofNoGlobalPetsc, isLocal);

          if (isLocal)
          {
            localDofsMapping[dofNoLocal] = iter->second;
          }
        }

        LOG(DEBUG) << "dofsMapping global: " << mapping.dofsMapping << " to local: " << localDofsMapping;

        // assign to dofs mapping
        mapping.dofsMapping = localDofsMapping;
      }

      // prepare communication if dofNoIsGlobalTo is global
      if (mapping.mode == DofsMappingType::modeCommunicate && mapping.dofNoIsGlobalTo)
      {
        std::map<int,std::vector<int>> remoteDofNosGlobalNaturalAtRanks;

        for (std::map<int,std::vector<int>>::iterator iter = mapping.dofsMapping.begin(); iter != mapping.dofsMapping.end(); iter++)
        {
          // if the "from" dofNo is local
          if (iter->first != -1)
          {
            // loop over the "to" dofNos
            for (int dofNoToGlobalNatural : iter->second)
            {
              // get the rank on which the dofNoToGlobalNatural is located
              int rankNo = meshPartitionBaseTo->getRankOfDofNoGlobalNatural(dofNoToGlobalNatural);

              // store the global dof no at the foreign rank
              remoteDofNosGlobalNaturalAtRanks[rankNo].push_back(dofNoToGlobalNatural);

              // store the local dof at the local rank
              mapping.dofNosLocalOfValuesToSendToRanks[rankNo].push_back(iter->first);
            }
          }
        }

        // initialize the communication data, send the global dof nos to the ranks where the data will be received
        std::vector<int> ownDofNosGlobalNatural;
        mapping.valueCommunicator.initialize(remoteDofNosGlobalNaturalAtRanks, ownDofNosGlobalNatural, data_.functionSpace()->meshPartition()->rankSubset());

        // convert the own received dof nos from global natural ordering to local numbering
        mapping.allDofNosToSetLocal.clear();

        for (int ownDofNoGlobalNatural : ownDofNosGlobalNatural)
        {
          global_no_t nodeNoGlobalNatural = ownDofNoGlobalNatural;  // this is not implemented for Hermite, we assume node=dof

          bool isOnLocalDomain;
          node_no_t nodeNoLocal = meshPartitionBaseTo->getNodeNoLocalFromGlobalNatural(nodeNoGlobalNatural, isOnLocalDomain);

          if (!isOnLocalDomain)
            LOG(FATAL) << "In MapDofs for mapping between connector slots " << mapping.connectorSlotNoFrom
              << " -> " << mapping.connectorSlotNosTo
              << ", node in global natural numbering " << nodeNoGlobalNatural << " is not on local domain.";

          dof_no_t dofNoLocal = nodeNoLocal;
          mapping.allDofNosToSetLocal.push_back(dofNoLocal);
        }
      }
      else
      {
        // no communication will take place
        // process dofsMapping
        for (std::map<int,std::vector<int>>::iterator iter = mapping.dofsMapping.begin(); iter != mapping.dofsMapping.end(); iter++)
        {
          // if the "from" dofNo is local
          if (iter->first != -1)
          {
            // loop over the "to" dofNos
            for (int dofNoTo : iter->second)
            {
              LOG(DEBUG) << "copy local dof " << iter->first << " -> " << dofNoTo;

              // if the dofs to map to are given as global natural no.s
              if (mapping.dofNoIsGlobalTo)
              {
                global_no_t nodeNoGlobalNatural = dofNoTo;  // this is not implemented for Hermite, we assume node=dof

                bool isOnLocalDomain;
                node_no_t nodeNoLocal = meshPartitionBaseTo->getNodeNoLocalFromGlobalNatural(nodeNoGlobalNatural, isOnLocalDomain);

                if (isOnLocalDomain)
                {
                  dof_no_t dofNoToLocal = nodeNoLocal;

                  // store the local dof at the local rank
                  LOG(DEBUG) << "dofNoTo was global, add dofNoFrom local " << iter->first << " for rank " << ownRankNo << ", now all: " 
                    << mapping.dofNosLocalOfValuesToSendToRanks;
                  mapping.dofNosLocalOfValuesToSendToRanks[ownRankNo].push_back(iter->first);
                  mapping.allDofNosToSetLocal.push_back(dofNoToLocal);
                }
              }
              else
              {
                // if the dofs to map to are given as local numbers
                dof_no_t dofNoToLocal = dofNoTo;

                // store the local dof at the local rank
                LOG(DEBUG) << "dofNoTo is local, add dofNoFrom local " << iter->first << " for rank " << ownRankNo << ", now all: " << mapping.dofNosLocalOfValuesToSendToRanks;
                mapping.dofNosLocalOfValuesToSendToRanks[ownRankNo].push_back(iter->first);
                mapping.allDofNosToSetLocal.push_back(dofNoToLocal);
              }
            }
          }
        }
      }

      // here, the following variables are set:
      // std::map<int,std::vector<dof_no_t>> dofNosLocalOfValuesToSendToRanks;      //< for every rank the local dof nos of the values that will be sent to the rank
      // std::vector<dof_no_t> allDofNosToSetLocal;                            //< for the received values the local dof nos where to store the values in the field variable

      LOG(DEBUG) << "initialized mapping slots " << mapping.connectorSlotNoFrom << " -> " << mapping.connectorSlotNosTo;
      LOG(DEBUG) << "  dofsMapping (local: " << (mapping.dofNoIsGlobalTo? "global" : "local") << ": " << mapping.dofsMapping;
      LOG(DEBUG) << "  dofNosLocalOfValuesToSendToRanks: " << mapping.dofNosLocalOfValuesToSendToRanks;
      LOG(DEBUG) << "  allDofNosToSetLocal:         " << mapping.allDofNosToSetLocal;
    }
  }  // loop over mappings
}

}  // namespace Control
