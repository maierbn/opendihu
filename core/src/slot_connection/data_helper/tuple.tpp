#include "slot_connection/data_helper/slot_connector_data_helper.h"

// ----------------------------------
// tuple

template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
int SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
nSlots(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData)
{
  return SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData))
    + SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));
}

//! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
int SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
nArrayItems(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData, int slotNo)
{
  if (!slotConnectorData)
    return 0;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    return SlotConnectorDataHelper<SlotConnectorDataType1>::nArrayItems(std::get<0>(*slotConnectorData), slotNo);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    return SlotConnectorDataHelper<SlotConnectorDataType2>::nArrayItems(std::get<1>(*slotConnectorData), offsetSlotNo);
  }
  return 0;
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

//! get the mesh name of the field variable at the slot
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
std::string SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
getMeshName(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo
)
{
  if (!slotConnectorData)
    return std::string();

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    return SlotConnectorDataHelper<SlotConnectorDataType1>::getMeshName(std::get<0>(*slotConnectorData), slotNo);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    return SlotConnectorDataHelper<SlotConnectorDataType2>::getMeshName(std::get<1>(*slotConnectorData), offsetSlotNo);
  }

  return std::string();
}

//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
bool SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
slotSetValues(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!slotConnectorData)
    return false;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    return SlotConnectorDataHelper<SlotConnectorDataType1>::slotSetValues(std::get<0>(*slotConnectorData), slotNo, arrayIndex, dofNosLocal, values, petscInsertMode);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    return SlotConnectorDataHelper<SlotConnectorDataType2>::slotSetValues(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values, petscInsertMode);
  }
  return true;
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

//! set the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
slotSetGeometryValues(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
    return;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    SlotConnectorDataHelper<SlotConnectorDataType1>::slotSetGeometryValues(std::get<0>(*slotConnectorData), slotNo, arrayIndex, dofNosLocal, values);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    SlotConnectorDataHelper<SlotConnectorDataType2>::slotSetGeometryValues(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values);
  }
}

//! get the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>::
slotGetGeometryValues(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
    return;

  int nSlotsFirstTuple = SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*slotConnectorData));
  int nSlotsSecondTuple = SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<1>(*slotConnectorData));

  if (slotNo < nSlotsFirstTuple)
  {
    SlotConnectorDataHelper<SlotConnectorDataType1>::slotGetGeometryValues(std::get<0>(*slotConnectorData), slotNo, arrayIndex, dofNosLocal, values);
  }
  else if (slotNo < nSlotsFirstTuple + nSlotsSecondTuple)
  {
    int offsetSlotNo = slotNo - nSlotsFirstTuple;
    SlotConnectorDataHelper<SlotConnectorDataType2>::slotGetGeometryValues(std::get<1>(*slotConnectorData), offsetSlotNo, arrayIndex, dofNosLocal, values);
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
