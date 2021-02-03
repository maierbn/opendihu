#include "slot_connection/data_helper/slot_connector_data_helper.h"

// ----------------------------------
// vector of vector

template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
nSlots(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData)
{
  if (!slotConnectorData)
    return 0;
  if (!slotConnectorData->empty())
    if (!(*slotConnectorData)[0]->empty())
      return (*(*slotConnectorData)[0])[0]->nSlots();
  return 0;
}

//! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
nArrayItems(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData)
{
  if (!slotConnectorData)
    return 0;
  if (!slotConnectorData->empty())
    if (!(*slotConnectorData)[0]->empty())
      return slotConnectorData->size() * (*slotConnectorData)[0]->size();
  return 0;
}

//! get the mesh partition of the field variable at the slot
template<typename SlotConnectorDataType>
std::shared_ptr<Partition::MeshPartitionBase> SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
getMeshPartitionBase(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return nullptr;
  if (slotConnectorData->empty())
    return nullptr;

  int sizeFirstVector = slotConnectorData->size();
  int sizeSecondVector = (*slotConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return nullptr;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*slotConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable2[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
}

//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
slotSetValues(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->empty())
    return;

  int sizeFirstVector = slotConnectorData->size();
  int sizeSecondVector = (*slotConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*slotConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].values;

    int componentNo = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << " in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable2[index].values;

    int componentNo = (*(*slotConnectorData)[arrayIndex1])->variable2[index].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
slotGetValues(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->empty())
    return;

  int sizeFirstVector = slotConnectorData->size();
  int sizeSecondVector = (*slotConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*slotConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].values;

    int componentNo = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable2[index].values;

    int componentNo = (*(*slotConnectorData)[arrayIndex1])->variable2[index].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! set the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
slotSetGeometryValues(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->empty())
    return;

  int sizeFirstVector = slotConnectorData->size();
  int sizeSecondVector = (*slotConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*slotConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].values;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "in fieldVariable \"" << fieldVariable->name() << "\", function space \"" << fieldVariable->functionSpace()->meshName() << "\""
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->functionSpace()->geometryField().setValues(dofNosLocal, values);
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable2[index].values;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "in fieldVariable \"" << fieldVariable->name() << "\", function space \"" << fieldVariable->functionSpace()->meshName() << "\""
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->functionSpace()->geometryField().setValues(dofNosLocal, values);
  }
}

//! get the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
slotGetGeometryValues(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->empty())
    return;

  int sizeFirstVector = slotConnectorData->size();
  int sizeSecondVector = (*slotConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*slotConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].values;

    // get the actual values
    fieldVariable->functionSpace()->geometryField().getValues(dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "from fieldVariable \"" << fieldVariable->name() << "\", function space \"" << fieldVariable->functionSpace()->meshName() << "\""
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable2[index].values;

    // get the actual values
    fieldVariable->functionSpace()->geometryField().getValues(dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "from fieldVariable \"" << fieldVariable->name() << "\", function space \"" << fieldVariable->functionSpace()->meshName() << "\""
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
slotSetRepresentationGlobal(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->empty())
    return;

  int sizeFirstVector = slotConnectorData->size();
  int sizeSecondVector = (*slotConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*slotConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable1[slotNo].values;

    fieldVariable->setRepresentationGlobal();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*slotConnectorData)[arrayIndex1])->variable2[index].values;

    fieldVariable->setRepresentationGlobal();
  }
}

//! collect all slot names
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
getSlotNames(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
  std::vector<std::string> &slotNames
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->empty())
    return;

  if (!(*slotConnectorData)[0]->empty())
  {
    // only collect the slot names from the first vector item, because all items should have the same slots
    SlotConnectorDataHelper<SlotConnectorDataType>::
      getSlotNames((*(*slotConnectorData)[0])[0], slotNames);
  }
}

//! get a string representation of the slot connector data for debugging
template<typename SlotConnectorDataType>
std::string SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>::
getString(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData)
{
  if (!slotConnectorData)
    return std::string("[[x]]");
  if (!slotConnectorData->empty())
    if (!(*slotConnectorData)[0]->empty())
    {
      std::stringstream s;
      s << "[";
      for (int j = 0; j < slotConnectorData->size(); j++)
      {
        if (j != 0)
          s << ", ";
        s << "[";
        for (int i = 0; i < (*slotConnectorData)[0]->size(); i++)
        {
          if (i != 0)
            s << ", ";
          s << SlotConnectorDataHelper<SlotConnectorDataType>::
            getString((*(*slotConnectorData)[0])[0]);
        }
        s << "]";
      }
      s << "]";

      return s.str();
    }
  return std::string("[[]]");
}
