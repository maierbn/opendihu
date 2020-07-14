#include "slot_connection/slot_connector_data_helper.h"

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

//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<SlotConnectorDataType>::
slotSetValues(
  std::shared_ptr<SlotConnectorDataType> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
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
  s << slotConnectorData;
  return s.str();
}

// vector

template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
nSlots(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData)
{
  if (!slotConnectorData)
    return 0;
  if (!slotConnectorData->empty())
    return (*slotConnectorData)[0]->nSlots();
  return 0;
}

//! get the mesh partition of the field variable at the slot
template<typename SlotConnectorDataType>
std::shared_ptr<Partition::MeshPartitionBase> SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
getMeshPartitionBase(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return nullptr;
  if (slotConnectorData->size() <= arrayIndex)
    return nullptr;

  int nSlotsVariable1 = (*slotConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable2[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
}


//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
slotSetValues(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->size() <= arrayIndex)
    return;

  int nSlotsVariable1 = (*slotConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable1[slotNo].values;

    int componentNo = (*slotConnectorData)[arrayIndex]->variable1[slotNo].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << ": "
      << "in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable2[index].values;

    int componentNo = (*slotConnectorData)[arrayIndex]->variable2[index].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << ": "
      << "in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
slotGetValues(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->size() <= arrayIndex)
    return;

  int nSlotsVariable1 = (*slotConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable1[slotNo].values;

    int componentNo = (*slotConnectorData)[arrayIndex]->variable1[slotNo].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << ": "
      << "from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable2[index].values;

    int componentNo = (*slotConnectorData)[arrayIndex]->variable2[index].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << ": "
      << "from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
slotSetRepresentationGlobal(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return;
  if (slotConnectorData->size() <= arrayIndex)
    return;

  int nSlotsVariable1 = (*slotConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable1Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable1[slotNo].values;

    fieldVariable->setRepresentationGlobal();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename SlotConnectorDataType::FieldVariable2Type> fieldVariable
      = (*slotConnectorData)[arrayIndex]->variable2[index].values;

    fieldVariable->setRepresentationGlobal();
  }
}

//! collect all slot names
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
getSlotNames(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  std::vector<std::string> &slotNames
)
{
  if (!slotConnectorData)
    return;

  if (!slotConnectorData->empty())
  {
    for (int i = 0; i < slotConnectorData->size(); i++)
    {
      SlotConnectorDataHelper<SlotConnectorDataType>::
        getSlotNames((*slotConnectorData)[i], slotNames);
    }
  }
}

//! get a string representation of the slot connector data for debugging
template<typename SlotConnectorDataType>
std::string SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
getString(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData)
{
  if (!slotConnectorData)
    return std::string("[x]");

  if (!slotConnectorData->empty())
  {
    std::stringstream s;
    s << "[";
    for (int i = 0; i < slotConnectorData->size(); i++)
    {
      if (i != 0)
        s << ", ";
      s << SlotConnectorDataHelper<SlotConnectorDataType>::
        getString((*slotConnectorData)[i]);
    }
    s << "]";
    return s.str();
  }
  return std::string("[]");
}

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
    for (int j = 0; j < slotConnectorData->size(); j++)
    {
      for (int i = 0; i < (*slotConnectorData)[0]->size(); i++)
      {
        SlotConnectorDataHelper<SlotConnectorDataType>::
          getSlotNames((*(*slotConnectorData)[0])[0], slotNames);
      }
    }
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

// tuple

template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
int SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
nSlots(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData)
{
  return SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData))
    + SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));
}

//! get the mesh partition of the field variable at the slot
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
std::shared_ptr<Partition::MeshPartitionBase> SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
getMeshPartitionBase(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return nullptr;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    return SlotConnectorDataHelper<SlotConnectorDataType1>::getMeshPartitionBase(std::get<0>(*slotConnectorData), slotNo, arrayIndex);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    return SlotConnectorDataHelper<SlotConnectorDataType2>::getMeshPartitionBase(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex);
  }

  return nullptr;
}

//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
slotSetValues(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!slotConnectorData)
    return;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    SlotConnectorDataHelper<SlotConnectorDataType1>::slotSetValues(std::get<0>(*slotConnectorData), slotNo, arrayIndex, dofNosLocal, values, petscInsertMode);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    SlotConnectorDataHelper<SlotConnectorDataType2>::slotSetValues(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values, petscInsertMode);
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
slotGetValues(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!slotConnectorData)
    return;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    SlotConnectorDataHelper<SlotConnectorDataType1>::slotGetValues(std::get<0>(*slotConnectorData), slotNo, arrayIndex, dofNosLocal, values);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    SlotConnectorDataHelper<SlotConnectorDataType2>::slotGetValues(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values);
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
slotSetRepresentationGlobal(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!slotConnectorData)
    return;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    SlotConnectorDataHelper<SlotConnectorDataType1>::slotSetRepresentationGlobal(std::get<0>(*slotConnectorData), slotNo, arrayIndex);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    SlotConnectorDataHelper<SlotConnectorDataType2>::slotSetRepresentationGlobal(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex);
  }
}

//! collect all slot names
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
getSlotNames(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  std::vector<std::string> &slotNames
)
{
  SlotConnectorDataHelper<SlotConnectorDataType1>::getSlotNames(std::get<0>(*slotConnectorData), slotNames);
  SlotConnectorDataHelper<SlotConnectorDataType2>::getSlotNames(std::get<1>(*slotConnectorData), slotNames);
}

//! get a string representation of the slot connector data for debugging
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
std::string SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
getString(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData)
{
  std::stringstream s;
  s << "{\n\t" << SlotConnectorDataHelper<SlotConnectorDataType1>::getString(std::get<0>(*slotConnectorData)) << "\n;\n\t"
    << SlotConnectorDataHelper<SlotConnectorDataType2>::getString(std::get<1>(*slotConnectorData)) << "\n}";
  return s.str();
}
