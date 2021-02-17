#pragma once

#include <Python.h>  // has to be the first included header

#include "partition/mesh_partition/00_mesh_partition_base.h"

/** Helper class that gives the number of connector slots and allows to access the slots' field variables.
 *  `arrayIndex` refers to the item of an array, e.g., resulting from MultipleInstances (e.g., multiple fibers)
 */
template<typename SlotConnectorDataType>
struct SlotConnectorDataHelper
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<SlotConnectorDataType> slotConnectorData);

  //! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
  static int nArrayItems(std::shared_ptr<SlotConnectorDataType> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static bool slotSetValues(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! set the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotSetGeometryValues(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values);

  //! get the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotGetGeometryValues(
    std::shared_ptr<SlotConnectorDataType> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
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

  //! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
  static int nArrayItems(std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static bool slotSetValues(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! set the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotSetGeometryValues(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values);

  //! get the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotGetGeometryValues(
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
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

  //! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
  static int nArrayItems(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static bool slotSetValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! set the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotSetGeometryValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values);

  //! get the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotGetGeometryValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType>>>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
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

  //! get the number of items if the slotConnector is organized in an array, `arrayIndex` can then be chosen in [0,nArrayItems]
  static int nArrayItems(std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static bool slotSetValues(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );

  //! set the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotSetGeometryValues(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<Vec3> &values);

  //! get the geometry values at given dofs of the function space of the field variable given by slotNo
  static void slotGetGeometryValues(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<Vec3> &values
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

#include "slot_connection/data_helper/pointer.tpp"
#include "slot_connection/data_helper/tuple.tpp"
#include "slot_connection/data_helper/vector.tpp"
#include "slot_connection/data_helper/vector_of_vector.tpp"
