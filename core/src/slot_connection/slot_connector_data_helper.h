#pragma once

#include <Python.h>  // has to be the first included header

#include "partition/mesh_partition/00_mesh_partition_base.h"

/** Helper class that gives the number of connector slots
 */
template<typename SlotConnectorDataType>
struct SlotConnectorDataHelper
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<SlotConnectorDataType> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! call setRepresentationGlobal on the field variable
  static void slotSetRepresentationGlobal(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! collect all slot names
  static void getSlotNames(std::shared_ptr<SlotConnectorDataType> slotConnectorData, std::vector<std::string> &slotNames);

  //! get a string representation of the slot connector data for debugging
  static std::string getString(std::shared_ptr<SlotConnectorDataType> slotConnectorData);
};

template<typename SlotConnectorDataType>
struct SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType>>>
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! call setRepresentationGlobal on the field variable
  static void slotSetRepresentationGlobal(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! collect all slot names, slots in the vector only get included once
  //! a, b, [(c,d),(c,d),(c,d)], e -> [a, b, c, d, e] and not [a, b, c, d, c, d, c, d, e]
  static void getSlotNames(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    std::vector<std::string> &slotNames
  );

  //! get a string representation of the slot connector data for debugging
  static std::string getString(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData);
};

template<typename SlotConnectorDataType>
struct SlotConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>>
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! call setRepresentationGlobal on the field variable
  static void slotSetRepresentationGlobal(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! collect all slot names, slots in the vector only get included once
  //! a, b, [(c,d),(c,d),(c,d)], e -> [a, b, c, d, e] and not [a, b, c, d, c, d, c, d, e]
  static void getSlotNames(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    std::vector<std::string> &slotNames
  );

  //! get a string representation of the slot connector data for debugging
  static std::string getString(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData);
};

template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
struct SlotConnectorDataHelper<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>>
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! call setRepresentationGlobal on the field variable
  static void slotSetRepresentationGlobal(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! collect all slot names
  static void getSlotNames(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    std::vector<std::string> &slotNames
  );

  //! get a string representation of the slot connector data for debugging
  static std::string getString(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData);
};

#include "slot_connection/slot_connector_data_helper.tpp"
