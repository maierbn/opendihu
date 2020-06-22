#pragma once

#include <Python.h>  // has to be the first included header

#include "partition/mesh_partition/00_mesh_partition_base.h"

/** Helper class that gives the number of output connector slots
 */
template<typename OutputConnectorDataType>
struct OutputConnectorDataHelper
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<OutputConnectorDataType> outputConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<OutputConnectorDataType> outputConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<OutputConnectorDataType> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<OutputConnectorDataType> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );
};

template<typename OutputConnectorDataType>
struct OutputConnectorDataHelper<std::vector<std::shared_ptr<OutputConnectorDataType>>>
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );
};

template<typename OutputConnectorDataType>
struct OutputConnectorDataHelper<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>>
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType>>>>> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );
};

template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
struct OutputConnectorDataHelper<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>>
{
  //! get the number of slots
  static int nSlots(std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData);

  //! get the mesh partition of the field variable at the slot
  static std::shared_ptr<Partition::MeshPartitionBase> getMeshPartitionBase(
    std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
    int slotNo, int arrayIndex
  );

  //! set the values at given dofs at the field variable given by slotNo
  static void slotSetValues(
    std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode = INSERT_VALUES
  );

  //! get the values at given dofs at the field variable given by slotNo
  static void slotGetValues(
    std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
    int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values
  );
};

#include "output_connector_data_transfer/output_connector_data_helper.tpp"
