#include "slot_connection/data_helper/slot_connector_data_helper.h"

// ----------------------------------
// vector

template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
nSlots(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData)
{
  if (!slotConnectorData)
    return 0;
  if (!slotConnectorData->empty())
    return SlotConnectorDataHelper<SlotConnectorDataType>::nSlots((*slotConnectorData)[0]);
  return 0;
}

//! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
template<typename SlotConnectorDataType>
int SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
nArrayItems(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData, int slotNo)
{
  if (!slotConnectorData)
    return 0;
  if (slotConnectorData->empty())
    return 0;

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  return slotConnectorData->size() * nArrayItemsPerItem;
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
  if (slotConnectorData->empty())
    return nullptr;

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  int itemNo = arrayIndex / nArrayItemsPerItem;
  int arrayIndexInItem = arrayIndex % nArrayItemsPerItem;

  return SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase((*slotConnectorData)[itemNo], slotNo, arrayIndexInItem);
}

//! set the values at given dofs at the field variable given by slotNo
template<typename SlotConnectorDataType>
bool SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
slotSetValues(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode
)
{
  if (!slotConnectorData)
  {
    LOG(DEBUG) << "slotSetValues (vector), slotConnectorData is not set";
    return false;
  }
  if (slotConnectorData->empty())
  {
    LOG(DEBUG) << "slotSetValues (vector), slotConnectorData is empty";
    return false;
  }

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  int itemNo = arrayIndex / nArrayItemsPerItem;
  int arrayIndexInItem = arrayIndex % nArrayItemsPerItem;

  return SlotConnectorDataHelper<SlotConnectorDataType>::slotSetValues((*slotConnectorData)[itemNo], slotNo, arrayIndexInItem,
                                                                       dofNosLocal, values, petscInsertMode);
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
  {
    LOG(DEBUG) << "slotGetValues (vector), slotConnectorData is not set";
    return;
  }
  if (slotConnectorData->empty())
  {
    LOG(DEBUG) << "slotGetValues (vector), slotConnectorData is empty";
    return;
  }

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  int itemNo = arrayIndex / nArrayItemsPerItem;
  int arrayIndexInItem = arrayIndex % nArrayItemsPerItem;

  SlotConnectorDataHelper<SlotConnectorDataType>::slotGetValues((*slotConnectorData)[itemNo], slotNo, arrayIndexInItem,
                                                                dofNosLocal, values);
}

//! set the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
slotSetGeometryValues(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
  {
    LOG(DEBUG) << "slotSetGeometryValues (vector), slotConnectorData is not set";
    return;
  }
  if (slotConnectorData->empty())
  {
    LOG(DEBUG) << "slotSetGeometryValues (vector), slotConnectorData is empty";
    return;
  }

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  int itemNo = arrayIndex / nArrayItemsPerItem;
  int arrayIndexInItem = arrayIndex % nArrayItemsPerItem;

  SlotConnectorDataHelper<SlotConnectorDataType>::slotSetGeometryValues((*slotConnectorData)[itemNo], slotNo, arrayIndexInItem,
                                                                        dofNosLocal, values);
}

//! get the geometry values at given dofs of the function space of the field variable given by slotNo
template<typename SlotConnectorDataType>
void SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>::
slotGetGeometryValues(
  std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
  int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
)
{
  if (!slotConnectorData)
  {
    LOG(DEBUG) << "slotGetGeometryValues (vector), slotConnectorData is not set";
    return;
  }
  if (slotConnectorData->empty())
  {
    LOG(DEBUG) << "slotGetGeometryValues (vector), slotConnectorData is empty";
    return;
  }

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  int itemNo = arrayIndex / nArrayItemsPerItem;
  int arrayIndexInItem = arrayIndex % nArrayItemsPerItem;

  LOG(DEBUG) << "slotGetGeometryValues on item " << arrayIndex << " -> " << arrayIndexInItem << "(" << nArrayItemsPerItem << ")";
  SlotConnectorDataHelper<SlotConnectorDataType>::slotGetGeometryValues((*slotConnectorData)[itemNo], slotNo, arrayIndexInItem,
                                                                        dofNosLocal, values);
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
  if (slotConnectorData->empty())
    return;

  int nArrayItemsPerItem = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems((*slotConnectorData)[0], slotNo);
  int itemNo = arrayIndex / nArrayItemsPerItem;
  int arrayIndexInItem = arrayIndex % nArrayItemsPerItem;

  SlotConnectorDataHelper<SlotConnectorDataType>::slotSetRepresentationGlobal((*slotConnectorData)[itemNo], slotNo, arrayIndexInItem);
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
    // only collect the slot names from the first item of the vector, as they should be all the same
    SlotConnectorDataHelper<SlotConnectorDataType>::
      getSlotNames((*slotConnectorData)[0], slotNames);
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
