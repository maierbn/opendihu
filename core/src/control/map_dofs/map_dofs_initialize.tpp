#include "control/map_dofs/map_dofs.h"

#include "output_connector_data_transfer/output_connector_data_helper.h"

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
  // parse settings
  int nAdditionalFieldVariables = this->specificSettings_.getOptionInt("nAdditionalFieldVariables", 0, PythonUtility::NonNegative);

  parseMappingFromSettings("beforeComputation", mappingsBeforeComputation_);
  parseMappingFromSettings("afterComputation", mappingsAfterComputation_);

  // initialize solver
  nestedSolver_.initialize();

  // create function space to use for the additional field variables
  std::shared_ptr<FunctionSpaceType> functionSpace = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);
  data_.setFunctionSpace(functionSpace);

  // initialize data, i.e. create additional field variables
  data_.initialize(nAdditionalFieldVariables, nestedSolver_.getOutputConnectorData());

  // prepare communication by initializing mappings
  initializeCommunication(mappingsBeforeComputation_);
  initializeCommunication(mappingsAfterComputation_);
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
parseMappingFromSettings(std::string settingsKey, std::vector<DofsMappingType> &mappings)
{
  PyObject *listPy = this->specificSettings_.getOptionPyObject("beforeComputation");
  std::vector<PyObject *> list = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);

  // loop over items of the list under "beforeComputation" or "afterComputation"
  for (int i = 0; i < list.size(); i++)
  {
    PythonConfig mappingSpecification(this->specificSettings_, i);

    DofsMappingType newDofsMapping;

    // parse options of this item
    newDofsMapping.outputConnectorSlotNoFrom = mappingSpecification.getOptionInt("fromOutputConnectorSlotNo", 0, PythonUtility::NonNegative);
    newDofsMapping.outputConnectorSlotNoTo   = mappingSpecification.getOptionInt("toOutputConnectorSlotNo",   0, PythonUtility::NonNegative);

    newDofsMapping.outputConnectorArrayIndexFrom = mappingSpecification.getOptionInt("fromOutputConnectorArrayIndex", 0, PythonUtility::NonNegative);
    newDofsMapping.outputConnectorArrayIndexTo   = mappingSpecification.getOptionInt("toOutputConnectorArrayIndex",   0, PythonUtility::NonNegative);

    newDofsMapping.dofNoIsGlobalFrom = mappingSpecification.getOptionBool("fromDofNosIsGlobal", true);
    newDofsMapping.dofNoIsGlobalTo   = mappingSpecification.getOptionBool("toDofNosIsGlobal",   false);

    std::string mode = mappingSpecification.getOptionString("mode", "copyLocal");
    if (mode == "copyLocal")
    {
      newDofsMapping.mode = DofsMappingType::modeCopyLocal;
    }
    else if (mode == "communicate")
    {
      newDofsMapping.mode = DofsMappingType::modeCommunicate;
    }
    else
    {
      LOG(FATAL) << mappingSpecification << "[\"mode\"] is \"" << mode << "\", but allowed values are \"copyLocal\" and \"communicate\".";
    }

    PyObject *object = mappingSpecification.getOptionPyObject("dofsMapping");
    newDofsMapping.dofsMapping = PythonUtility::convertFromPython<std::map<int,std::vector<int>>>::get(object);
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
std::shared_ptr<Partition::MeshPartitionBase> MapDofs<FunctionSpaceType,NestedSolverType>::
getMeshPartitionBase(int slotNo, int arrayIndex)
{
  /*
  typedef std::tuple<
    std::shared_ptr<typename NestedSolverType::OutputConnectorDataType>,
    std::shared_ptr<OutputConnectorData<FunctionSpaceType,1>>
  > OutputConnectorDataType;
  */

  // need function space of affected field variables
  int nSlotsNestedSolver = OutputConnectorDataHelper<typename NestedSolverType::OutputConnectorDataType>::nSlots(
    std::get<0>(*data_.getOutputConnectorData())
  );

  if (slotNo < nSlotsNestedSolver)
  {
    return OutputConnectorDataHelper<typename NestedSolverType::OutputConnectorDataType>::getMeshPartitionBase(
      std::get<0>(*data_.getOutputConnectorData()), slotNo, arrayIndex
    );
  }
  else
  {
    int index = slotNo - nSlotsNestedSolver;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
      = std::get<1>(*data_.getOutputConnectorData())->variable1[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }

  return nullptr;
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
slotGetValues(int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values)
{
  // need function space of affected field variables
  int nSlotsNestedSolver = OutputConnectorDataHelper<typename NestedSolverType::OutputConnectorDataType>::nSlots(
    std::get<0>(*data_.getOutputConnectorData())
  );
  if (slotNo < nSlotsNestedSolver)
  {
    OutputConnectorDataHelper<typename NestedSolverType::OutputConnectorDataType>::slotGetValues(
      std::get<0>(*data_.getOutputConnectorData()), slotNo, arrayIndex, dofNosLocal, values
    );
  }
  else
  {
    // get the field variable and component no for this slot
    int index = slotNo - nSlotsNestedSolver;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
      = std::get<1>(*data_.getOutputConnectorData())->variable1[index].values;

    int componentNo = std::get<1>(*data_.getOutputConnectorData())->variable1[index].componentNo;  // should be 0

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
slotSetValues(int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode)
{
  // need function space of affected field variables
  int nSlotsNestedSolver = OutputConnectorDataHelper<typename NestedSolverType::OutputConnectorDataType>::nSlots(
    std::get<0>(*data_.getOutputConnectorData())
  );
  if (slotNo < nSlotsNestedSolver)
  {
    OutputConnectorDataHelper<typename NestedSolverType::OutputConnectorDataType>::slotSetValues(
      std::get<0>(*data_.getOutputConnectorData()), slotNo, arrayIndex, dofNosLocal, values, petscInsertMode
    );
  }
  else
  {
    // get the field variable and component no for this slot
    int index = slotNo - nSlotsNestedSolver;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> fieldVariable
      = std::get<1>(*data_.getOutputConnectorData())->variable1[index].values;

    int componentNo = std::get<1>(*data_.getOutputConnectorData())->variable1[index].componentNo;  // should be 0

    LOG(DEBUG) << "slot " << slotNo << ": in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
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
    std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBaseFrom = getMeshPartitionBase(mapping.outputConnectorSlotNoFrom, mapping.outputConnectorArrayIndexFrom);

    // get the mesh partition of the "to" field variable
    std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBaseTo = getMeshPartitionBase(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo);

    if (!meshPartitionBaseFrom)
      LOG(FATAL) << "MapDofs: Could not get mesh partition for \"from\" function space for mapping " << mapping.outputConnectorSlotNoFrom << " -> " << mapping.outputConnectorSlotNoTo;
    if (!meshPartitionBaseTo)
      LOG(FATAL) << "MapDofs: Could not get mesh partition for \"to\" function space for mapping " << mapping.outputConnectorSlotNoFrom << " -> " << mapping.outputConnectorSlotNoTo;

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
      mapping.receivedValueDofNosLocal.clear();

      for (int ownDofNoGlobalNatural : ownDofNosGlobalNatural)
      {
        global_no_t nodeNoGlobalNatural = ownDofNoGlobalNatural;  // this is not implemented for Hermite, we assume node=dof

        bool isOnLocalDomain;
        node_no_t nodeNoLocal = meshPartitionBaseTo->getNodeNoLocalFromGlobalNatural(nodeNoGlobalNatural, isOnLocalDomain);

        if (!isOnLocalDomain)
          LOG(FATAL) << "In MapDofs for mapping between output connector slots " << mapping.outputConnectorSlotNoFrom
            << " -> " << mapping.outputConnectorSlotNoTo
            << ", node in global natural numbering " << nodeNoGlobalNatural << " is not on local domain.";

        dof_no_t dofNoLocal = nodeNoLocal;
        mapping.receivedValueDofNosLocal.push_back(dofNoLocal);
      }
    }
    else
    {
      // no communication will take place
      for (std::map<int,std::vector<int>>::iterator iter = mapping.dofsMapping.begin(); iter != mapping.dofsMapping.end(); iter++)
      {
        // if the "from" dofNo is local
        if (iter->first != -1)
        {
          // loop over the "to" dofNos
          for (int dofNoTo : iter->second)
          {
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
                mapping.dofNosLocalOfValuesToSendToRanks[ownRankNo].push_back(iter->first);
                mapping.receivedValueDofNosLocal.push_back(dofNoToLocal);
              }
            }
            else
            {
              // if the dofs to map to are given as local numbers
              dof_no_t dofNoToLocal = dofNoTo;

              // store the local dof at the local rank
              mapping.dofNosLocalOfValuesToSendToRanks[ownRankNo].push_back(iter->first);
              mapping.receivedValueDofNosLocal.push_back(dofNoToLocal);
            }
          }
        }
      }
    }

    // here, the following variables are set:
    // std::map<int,std::vector<dof_no_t>> dofNosLocalOfValuesToSendToRanks;      //< for every rank the local dof nos of the values that will be sent to the rank
    // std::vector<dof_no_t> receivedValueDofNosLocal;                            //< for the received values the local dof nos where to store the values in the field variable

    LOG(DEBUG) << "initialized mapping slots " << mapping.outputConnectorSlotNoFrom << " -> " << mapping.outputConnectorSlotNoTo;
    LOG(DEBUG) << "  dofsMapping (local: " << (mapping.dofNoIsGlobalTo? "global" : "local") << ": " << mapping.dofsMapping;
    LOG(DEBUG) << "  dofNosLocalOfValuesToSendToRanks: " << mapping.dofNosLocalOfValuesToSendToRanks;
    LOG(DEBUG) << "  receivedValueDofNosLocal:         " << mapping.receivedValueDofNosLocal;

  }  // loop over mappings
}

}  // namespace Control
