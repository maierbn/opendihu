#include "slot_connection/data_helper/slot_connector_data_helper.h"

/** Helper class that gives the number of connector slots
 */
template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<SlotConnectorDataType>::
nSlots(std::shared_ptr<SlotConnectorDataType> slotConnectorData)
{
  if (!slotConnectorData)
    return 0;
  return slotConnectorData->nSlots();
}

//! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<SlotConnectorDataType>::
nArrayItems(std::shared_ptr<SlotConnectorDataType> slotConnectorData, int slotNo)
{
  if (!slotConnectorData)
    return 0;
  return 1;
}

//! get the mesh partition of the field variable at the slot
template<typename SlotConnectorDataType>
std::shared_ptr<Partition::MeshPartitionBase> SlotConnectorDataHelper<SlotConnectorDataType>::
getMeshPartitionBase(std::shared_ptr<SlotConnectorDataType> slotConnectorData, int slotNo, int arrayIndex)
{
  if (!slotConnectorData)
    return nullptr;
  int nSlotsVariable1 = slotConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
}

//! get the mesh name of the field variable at the slot
template<typename SlotConnectorDataType>
std::string SlotConnectorDataHelper<SlotConnectorDataType>::
getMeshName(std::shared_ptr<SlotConnectorDataType> slotConnectorData, int slotNo)
{
  if (!slotConnectorData)
    return std::string();
  int nSlotsVariable1 = slotConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshName();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;

    return fieldVariable->functionSpace()->meshName();
  }
  return std::string();
}

//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
bool SlotConnectorDataHelper<SlotConnectorDataType>::
slotSetValues(
  std::shared_ptr<SlotConnectorDataType> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!slotConnectorData)
    return false;

  int nSlotsVariable1 = slotConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    int componentNo = slotConnectorData->variable1[slotNo].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ": in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;

    int componentNo = slotConnectorData->variable2[index].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ": in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
  return true;
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<SlotConnectorDataType>::
slotGetValues(
  std::shared_ptr<SlotConnectorDataType> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!slotConnectorData)
    return;

  int nSlotsVariable1 = slotConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    int componentNo = slotConnectorData->variable1[slotNo].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;

    int componentNo = slotConnectorData->variable2[index].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! set the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<SlotConnectorDataType>::
slotSetGeometryValues(
  std::shared_ptr<SlotConnectorDataType> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
    return;

  int nSlotsVariable1 = slotConnectorData->variable1.size();
  using FunctionSpaceType = typename SlotConnectorDataType::FieldVariable1Type::FunctionSpace;
  using GeometryFieldType = FieldVariable::FieldVariable<FunctionSpaceType,3>;

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ": in fieldVariable \"" << fieldVariable->name() << "\", "
      << "function space \"" << fieldVariable->functionSpace()->meshName() << "\", "
      << "set dofs " << dofNosLocal << " of geometry field to values " << values;
    fieldVariable->functionSpace()->geometryField().setValues(dofNosLocal, values);
    
    // add the geometry field in the slot connector data, such that it will be automatically transferred to the connected slots
    slotConnectorData->addGeometryField(std::make_shared<GeometryFieldType>(fieldVariable->functionSpace()->geometryField()));
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ": in fieldVariable \"" << fieldVariable->name() << "\", "
      << "function space \"" << fieldVariable->functionSpace()->meshName() << "\", "
      << "set dofs " << dofNosLocal << " of geometry field to values " << values;
    fieldVariable->functionSpace()->geometryField().setValues(dofNosLocal, values);

    // add the geometry field in the slot connector data, such that it will be automatically transferred to the connected slots
    slotConnectorData->addGeometryField(std::make_shared<GeometryFieldType>(fieldVariable->functionSpace()->geometryField()));
  }
}

//! get the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<SlotConnectorDataType>::
slotGetGeometryValues(
  std::shared_ptr<SlotConnectorDataType> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
  {
    LOG(DEBUG) << "slotGetGeometryValues, slotConnectorData is not set";
    return;
  }

  int nSlotsVariable1 = slotConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    // get the actual values
    fieldVariable->functionSpace()->geometryField().getValues(dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from geometry field of fieldVariable \"" << fieldVariable->name() << "\", "
      << "function space \"" << fieldVariable->functionSpace()->meshName() << "\", "
      << "at dofs " << dofNosLocal << " get values " << values;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;

    // get the actual values
    fieldVariable->functionSpace()->geometryField().getValues(dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from geometry field of fieldVariable \"" << fieldVariable->name() << "\", "
      << "function space \"" << fieldVariable->functionSpace()->meshName() << "\", "
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<SlotConnectorDataType>::
slotSetRepresentationGlobal(
  std::shared_ptr<SlotConnectorDataType> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return;

  int nSlotsVariable1 = slotConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = slotConnectorData->variable1[slotNo].values;

    fieldVariable->setRepresentationGlobal();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = slotConnectorData->variable2[index].values;
  }
}

//! collect all slot names
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<SlotConnectorDataType>::
getSlotNames(std::shared_ptr<SlotConnectorDataType> slotConnectorData, std::vector<std::string> &slotNames)
{
  if (!slotConnectorData)
    return;

  slotNames.insert(slotNames.end(), slotConnectorData->slotNames.begin(), slotConnectorData->slotNames.end());
}

//! get a string representation of the slot connector data for debugging
template<typename SlotConnectorDataType>
std::string SlotConnectorDataHelper<SlotConnectorDataType>::
getString(std::shared_ptr<SlotConnectorDataType> slotConnectorData)
{
  if (!slotConnectorData)
    return std::string("x");

  std::stringstream s;
  s << *slotConnectorData;
  return s.str();
}
