#include "output_connector_data_transfer/output_connector_data_helper.h"

/** Helper class that gives the number of output connector slots
 */
template<typename OutputConnectorDataType>
int OutputConnectorDataHelper<OutputConnectorDataType>::
nSlots(std::shared_ptr<OutputConnectorDataType> outputConnectorData)
{
  if (!outputConnectorData)
    return 0;
  return outputConnectorData->nSlots();
}

//! get the mesh partition of the field variable at the slot
template<typename OutputConnectorDataType>
std::shared_ptr<Partition::MeshPartitionBase> OutputConnectorDataHelper<OutputConnectorDataType>::
getMeshPartitionBase(std::shared_ptr<OutputConnectorDataType> outputConnectorData, int slotNo, int arrayIndex)
{
  if (!outputConnectorData)
    return nullptr;
  int nSlotsVariable1 = outputConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = outputConnectorData->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = outputConnectorData->variable2[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
}

//! set the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<OutputConnectorDataType>::
slotSetValues(
  std::shared_ptr<OutputConnectorDataType> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!outputConnectorData)
    return;

  int nSlotsVariable1 = outputConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = outputConnectorData->variable1[slotNo].values;

    int componentNo = outputConnectorData->variable1[slotNo].componentNo;

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
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = outputConnectorData->variable2[index].values;

    int componentNo = outputConnectorData->variable2[index].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ": in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<OutputConnectorDataType>::
slotGetValues(
  std::shared_ptr<OutputConnectorDataType> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!outputConnectorData)
    return;

  int nSlotsVariable1 = outputConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = outputConnectorData->variable1[slotNo].values;

    int componentNo = outputConnectorData->variable1[slotNo].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = outputConnectorData->variable2[index].values;

    int componentNo = outputConnectorData->variable2[index].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<OutputConnectorDataType>::
slotSetRepresentationGlobal(
  std::shared_ptr<OutputConnectorDataType> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return;

  int nSlotsVariable1 = outputConnectorData->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = outputConnectorData->variable1[slotNo].values;

    fieldVariable->setRepresentationGlobal();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = outputConnectorData->variable2[index].values;
  }
}

//! collect all slot names
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<OutputConnectorDataType>::
getSlotNames(std::shared_ptr<OutputConnectorDataType> outputConnectorData, std::vector<std::string> &slotNames)
{
  if (!outputConnectorData)
    return;

  slotNames.insert(slotNames.end(), outputConnectorData->slotNames.begin(), outputConnectorData->slotNames.end());
}

//! get a string representation of the output connector data for debugging
template<typename OutputConnectorDataType>
std::string OutputConnectorDataHelper<OutputConnectorDataType>::
getString(std::shared_ptr<OutputConnectorDataType> outputConnectorData)
{
  if (!outputConnectorData)
    return std::string("x");

  std::stringstream s;
  s << outputConnectorData;
  return s.str();
}

// vector

template<typename OutputConnectorDataType>
int OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
nSlots(std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData)
{
  if (!outputConnectorData)
    return 0;
  if (!outputConnectorData->empty())
    return (*outputConnectorData)[0]->nSlots();
  return 0;
}

//! get the mesh partition of the field variable at the slot
template<typename OutputConnectorDataType>
std::shared_ptr<Partition::MeshPartitionBase> OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
getMeshPartitionBase(
  std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return nullptr;
  if (outputConnectorData->size() <= arrayIndex)
    return nullptr;

  int nSlotsVariable1 = (*outputConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable2[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
}


//! set the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
slotSetValues(
  std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->size() <= arrayIndex)
    return;

  int nSlotsVariable1 = (*outputConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable1[slotNo].values;

    int componentNo = (*outputConnectorData)[arrayIndex]->variable1[slotNo].componentNo;

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
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable2[index].values;

    int componentNo = (*outputConnectorData)[arrayIndex]->variable2[index].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << ": "
      << "in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
slotGetValues(
  std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->size() <= arrayIndex)
    return;

  int nSlotsVariable1 = (*outputConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable1[slotNo].values;

    int componentNo = (*outputConnectorData)[arrayIndex]->variable1[slotNo].componentNo;

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
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable2[index].values;

    int componentNo = (*outputConnectorData)[arrayIndex]->variable2[index].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << ": "
      << "from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
slotSetRepresentationGlobal(
  std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->size() <= arrayIndex)
    return;

  int nSlotsVariable1 = (*outputConnectorData)[arrayIndex]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable1[slotNo].values;

    fieldVariable->setRepresentationGlobal();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*outputConnectorData)[arrayIndex]->variable2[index].values;

    fieldVariable->setRepresentationGlobal();
  }
}

//! collect all slot names
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
getSlotNames(
  std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
  std::vector<std::string> &slotNames
)
{
  if (!outputConnectorData)
    return;

  if (!outputConnectorData->empty())
  {
    for (int i = 0; i < outputConnectorData->size(); i++)
    {
      OutputConnectorDataHelper<OutputConnectorDataType>::
        getSlotNames((*outputConnectorData)[i], slotNames);
    }
  }
}

//! get a string representation of the output connector data for debugging
template<typename OutputConnectorDataType>
std::string OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>::
getString(std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData)
{
  if (!outputConnectorData)
    return std::string("[x]");

  if (!outputConnectorData->empty())
  {
    std::stringstream s;
    s << "[";
    for (int i = 0; i < outputConnectorData->size(); i++)
    {
      if (i != 0)
        s << ", ";
      s << OutputConnectorDataHelper<OutputConnectorDataType>::
        getString((*outputConnectorData)[i]);
    }
    s << "]";
    return s.str();
  }
  return std::string("[]");
}

// vector of vector

template<typename OutputConnectorDataType>
int OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
nSlots(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData)
{
  if (!outputConnectorData)
    return 0;
  if (!outputConnectorData->empty())
    if (!(*outputConnectorData)[0]->empty())
      return (*(*outputConnectorData)[0])[0]->nSlots();
  return 0;
}

//! get the mesh partition of the field variable at the slot
template<typename OutputConnectorDataType>
std::shared_ptr<Partition::MeshPartitionBase> OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
getMeshPartitionBase(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return nullptr;
  if (outputConnectorData->empty())
    return nullptr;

  int sizeFirstVector = outputConnectorData->size();
  int sizeSecondVector = (*outputConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return nullptr;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*outputConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable1[slotNo].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable2[index].values;

    return fieldVariable->functionSpace()->meshPartitionBase();
  }
}

//! set the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
slotSetValues(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->empty())
    return;

  int sizeFirstVector = outputConnectorData->size();
  int sizeSecondVector = (*outputConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*outputConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable1[slotNo].values;

    int componentNo = (*(*outputConnectorData)[arrayIndex1])->variable1[slotNo].componentNo;

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
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable2[index].values;

    int componentNo = (*(*outputConnectorData)[arrayIndex1])->variable2[index].componentNo;

    // set the values in the field variable
    LOG(DEBUG) << "slot " << slotNo << ", array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "in fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", set dofs " << dofNosLocal << " to values " << values;
    fieldVariable->setValues(componentNo, dofNosLocal, values, petscInsertMode);
    LOG(DEBUG) << "fieldVariable: " << *fieldVariable;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
slotGetValues(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->empty())
    return;

  int sizeFirstVector = outputConnectorData->size();
  int sizeSecondVector = (*outputConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*outputConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable1[slotNo].values;

    int componentNo = (*(*outputConnectorData)[arrayIndex1])->variable1[slotNo].componentNo;

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
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable2[index].values;

    int componentNo = (*(*outputConnectorData)[arrayIndex1])->variable2[index].componentNo;

    // get the actual values
    fieldVariable->getValues(componentNo, dofNosLocal, values);

    LOG(DEBUG) << "slot " << slotNo << ": array index " << arrayIndex << " = [" << arrayIndex1 << "," << arrayIndex2 << "]: "
      << "from fieldVariable \"" << fieldVariable->name() << "\", component " << componentNo
      << ", at dofs " << dofNosLocal << " get values " << values;
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
slotSetRepresentationGlobal(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->empty())
    return;

  int sizeFirstVector = outputConnectorData->size();
  int sizeSecondVector = (*outputConnectorData)[0]->size();

  if (arrayIndex >= sizeFirstVector*sizeSecondVector)
    return;

  int arrayIndex1 = arrayIndex / sizeSecondVector;
  int arrayIndex2 = arrayIndex % sizeSecondVector;

  int nSlotsVariable1 = (*(*outputConnectorData)[arrayIndex1])[arrayIndex2]->variable1.size();

  // if the slot no corresponds to a field variables stored under variable1
  if (slotNo < nSlotsVariable1)
  {
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable1Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable1[slotNo].values;

    fieldVariable->setRepresentationGlobal();
  }
  else
  {
    // if the slot no corresponds to a field variables stored under variable2
    int index = slotNo - nSlotsVariable1;
    std::shared_ptr<typename OutputConnectorDataType::FieldVariable2Type> fieldVariable
      = (*(*outputConnectorData)[arrayIndex1])->variable2[index].values;

    fieldVariable->setRepresentationGlobal();
  }
}

//! collect all slot names
template<typename OutputConnectorDataType>
void OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
getSlotNames(
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
  std::vector<std::string> &slotNames
)
{
  if (!outputConnectorData)
    return;
  if (outputConnectorData->empty())
    return;

  if (!(*outputConnectorData)[0]->empty())
  {
    for (int j = 0; j < outputConnectorData->size(); j++)
    {
      for (int i = 0; i < (*outputConnectorData)[0]->size(); i++)
      {
        OutputConnectorDataHelper<OutputConnectorDataType>::
          getSlotNames((*(*outputConnectorData)[0])[0], slotNames);
      }
    }
  }
}

//! get a string representation of the output connector data for debugging
template<typename OutputConnectorDataType>
std::string OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>::
getString(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData)
{
  if (!outputConnectorData)
    return std::string("[[x]]");
  if (!outputConnectorData->empty())
    if (!(*outputConnectorData)[0]->empty())
    {
      std::stringstream s;
      s << "[";
      for (int j = 0; j < outputConnectorData->size(); j++)
      {
        if (j != 0)
          s << ", ";
        s << "[";
        for (int i = 0; i < (*outputConnectorData)[0]->size(); i++)
        {
          if (i != 0)
            s << ", ";
          s << OutputConnectorDataHelper<OutputConnectorDataType>::
            getString((*(*outputConnectorData)[0])[0]);
        }
        s << "]";
      }
      s << "]";

      return s.str();
    }
  return std::string("[[]]");
}

// tuple

template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
int OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
nSlots(std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData)
{
  return OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*outputConnectorData))
    + OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<1>(*outputConnectorData));
}

//! get the mesh partition of the field variable at the slot
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
std::shared_ptr<Partition::MeshPartitionBase> OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
getMeshPartitionBase(
  std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return nullptr;

  int nSlotsFirstTuple = OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*outputConnectorData));
  int nSlotsSecondTuple = OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<1>(*outputConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    return OutputConnectorDataHelper<OutputConnectorDataType1>::getMeshPartitionBase(std::get<0>(*outputConnectorData), slotNo, arrayIndex);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    return OutputConnectorDataHelper<OutputConnectorDataType2>::getMeshPartitionBase(std::get<1>(*outputConnectorData), offsetSlotNo, arrayIndex);
  }

  return nullptr;
}

//! set the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
slotSetValues(
  std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!outputConnectorData)
    return;

  int nSlotsFirstTuple = OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*outputConnectorData));
  int nSlotsSecondTuple = OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<1>(*outputConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    OutputConnectorDataHelper<OutputConnectorDataType1>::slotSetValues(std::get<0>(*outputConnectorData), slotNo, arrayIndex, dofNosLocal, values, petscInsertMode);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    OutputConnectorDataHelper<OutputConnectorDataType2>::slotSetValues(std::get<1>(*outputConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values, petscInsertMode);
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
slotGetValues(
  std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
)
{
  if (!outputConnectorData)
    return;

  int nSlotsFirstTuple = OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*outputConnectorData));
  int nSlotsSecondTuple = OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<1>(*outputConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    OutputConnectorDataHelper<OutputConnectorDataType1>::slotGetValues(std::get<0>(*outputConnectorData), slotNo, arrayIndex, dofNosLocal, values);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    OutputConnectorDataHelper<OutputConnectorDataType2>::slotGetValues(std::get<1>(*outputConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values);
  }
}

//! get the values at given dofs at the field variable given by slotNo
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
slotSetRepresentationGlobal(
  std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
  int slotNo, int arrayIndex
)
{
  if (!outputConnectorData)
    return;

  int nSlotsFirstTuple = OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*outputConnectorData));
  int nSlotsSecondTuple = OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<1>(*outputConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    OutputConnectorDataHelper<OutputConnectorDataType1>::slotSetRepresentationGlobal(std::get<0>(*outputConnectorData), slotNo, arrayIndex);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    OutputConnectorDataHelper<OutputConnectorDataType2>::slotSetRepresentationGlobal(std::get<1>(*outputConnectorData), offsetSlotNo, arrayIndex);
  }
}

//! collect all slot names
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
getSlotNames(
  std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
  std::vector<std::string> &slotNames
)
{
  OutputConnectorDataHelper<OutputConnectorDataType1>::getSlotNames(std::get<0>(*outputConnectorData), slotNames);
  OutputConnectorDataHelper<OutputConnectorDataType2>::getSlotNames(std::get<1>(*outputConnectorData), slotNames);
}

//! get a string representation of the output connector data for debugging
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
std::string OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>::
getString(std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData)
{
  std::stringstream s;
  s << "{\n\t" << OutputConnectorDataHelper<OutputConnectorDataType1>::getString(std::get<0>(*outputConnectorData)) << "\n;\n\t"
    << OutputConnectorDataHelper<OutputConnectorDataType2>::getString(std::get<1>(*outputConnectorData)) << "\n}";
  return s.str();
}
