#include "slot_connection/slot_connector_data_transfer_vector.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

#include "slot_connection/data_helper/slot_connector_data_helper.h"

/** Transfer between two vectors of any type
 */
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void SlotConnectorDataTransfer<
  std::vector<std::shared_ptr<SlotConnectorDataType1>>,
  std::vector<std::shared_ptr<SlotConnectorDataType2>>
>::transfer(const std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType1>>> transferableSolutionData1,
            std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorDataType2>>> transferableSolutionData2,
            SlotsConnection &slotsConnection,
            int offsetSlotNoData1, int offsetSlotNoData2,
            bool originalTransferDirection)
{
  // debugging and error output
  LOG(DEBUG) << "transfer vector (1)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "offsetSlot: " << offsetSlotNoData1 << ", " << offsetSlotNoData2 << ", slotsConnection: " << slotsConnection.getDebugInformation();
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType1>>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<std::vector<std::shared_ptr<SlotConnectorDataType2>>>::getString(transferableSolutionData2);
#endif



  VLOG(1) << "Solution vector mapping (slot_connector_data_transfer_vector.tpp)";
  int nTransferableVariables = std::min(transferableSolutionData1->size(), transferableSolutionData2->size());
  if (transferableSolutionData1->size() != transferableSolutionData2->size())
  {
    LOG(ERROR) << "Trying to transfer data from " << transferableSolutionData1->size() << " variables to " << transferableSolutionData2->size() << ", number has to be equal. "
      << "Now only using the first " << nTransferableVariables << " variable" << (nTransferableVariables != 1? "s" : "") << ". Types: " << std::endl
      << StringUtility::demangle(typeid(SlotConnectorDataType1).name()) << std::endl
      << StringUtility::demangle(typeid(SlotConnectorDataType2).name());
  }

  // loop over the two vectors A and B and transfer A[i] -> B[i] for every component
  for (int i = 0; i < nTransferableVariables; i++)
  {
    SlotConnectorDataTransfer<SlotConnectorDataType1,SlotConnectorDataType2>::transfer(
      (*transferableSolutionData1)[i],
      (*transferableSolutionData2)[i],
      slotsConnection,
      offsetSlotNoData1, offsetSlotNoData2,
      originalTransferDirection
    );
  }

  VLOG(1) << "at the end of std::vector transfer: transferableSolutionData2: " << (*transferableSolutionData2);
}

/** Transfer between a normal entry and a vector, the first entry of the vector is used
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void SlotConnectorDataTransfer<
  std::vector<std::shared_ptr<
    Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>
  >>,
  Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
>::transfer(const std::shared_ptr<std::vector<std::shared_ptr<
              Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>
            >>> transferableSolutionData1,
            std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> transferableSolutionData2,
            SlotsConnection &slotsConnection,
            int offsetSlotNoData1, int offsetSlotNoData2,
            bool originalTransferDirection)
{
  // debugging and error output
  LOG(DEBUG) << "transfer vector (2)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "offsetSlot: " << offsetSlotNoData1 << ", " << offsetSlotNoData2 << ", slotsConnection: " << slotsConnection.getDebugInformation();
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<std::vector<std::shared_ptr<
      Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>
    >>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>>::getString(transferableSolutionData2);
#endif

  VLOG(1) << "Solution vector mapping (slot_connector_data_transfer_vector.tpp)";

  if (transferableSolutionData1->empty())
  {
    LOG(ERROR) << "Trying to transfer data from an empty entry to a single vector. Types: " << std::endl
      << StringUtility::demangle(typeid(Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>).name()) << std::endl
      << StringUtility::demangle(typeid(Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>).name());
    return;
  }

  /** Explanation: map values from transferableSolutionData1 to transferableSolutionData2.
   *  These two variables specify all slots of two terms. They contain multiple field variables and the specification of component no's.
   *  "transferableSolutionData1" is of type (roughly) vector<field variables>. This results e.g. from multiple fibers.
   *  "transferableSolutionData2" is just normal field variables without a vector.
   *  Each object of type SlotConnectorData contains possibly two different types of field variables, variable1 and variable2.
   *  For cellML, these are for storing states and algebraics.
   *  Now this method maps components of field variables (i.e., scalar field variables) from 1 to 2 according to the slot connections
   *  given by "slotsConnection".
   *  We iterate over "1".variable1 and transfer data to "2".variable1 and "2".variable2. Then we iterate over "1".variable2
   *  and again map data to "2".variable1 and "2".variable2.
   */

  /*
  SlotConnectorDataTransfer<SlotConnectorDataType1,SlotConnectorDataType2>::transfer(
    (*transferableSolutionData1)[0],
    transferableSolutionData2,
    slotsConnection
  );
  */

  assert(transferableSolutionData1->size() > 0);

  std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1Front
    = (*transferableSolutionData1)[0];

  // initialize output connection object
  slotsConnection.initialize(*transferableSolutionData1Front, *transferableSolutionData2, offsetSlotNoData1, offsetSlotNoData2);

  // for the first vector of variables (the "states" in case of CellMLAdapter)
  for (int i = 0; i < transferableSolutionData1Front->variable1.size(); i++)
  {
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = slotsConnection.getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a> FieldVariable1;

    std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1Front->variable1[fromVectorIndex].values;
    int componentNo1                               = transferableSolutionData1Front->variable1[fromVectorIndex].componentNo;

    LOG(DEBUG) << "map slot from variable1, index " << fromVectorIndex << " (" << fieldVariable1->name() << "[" << componentNo1 << "])"
      << " to variable" << toVectorNo+1 << ", index " << toVectorIndex;

    if (componentNo1 < 0)
    {
      LOG(DEBUG) << "do not map this slot";
      continue;
    }

    if (toVectorNo == 0)
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable1[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2->variable1[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "] (" << fieldVariable2 << "), avoidCopyIfPossible: " << avoidCopyIfPossible << "(5)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData1->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable1 to transferableSolutionData2->variable1

        typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a> FieldVariable1;
        std::shared_ptr<FieldVariable1> fieldVariable1 = (*transferableSolutionData1)[fiberIndex]->variable1[fromVectorIndex].values;
        int componentNo1                               = (*transferableSolutionData1)[fiberIndex]->variable1[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
    else
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable2[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2->variable2[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(6)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData1->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable1 to transferableSolutionData2->variable2

        typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a> FieldVariable1;
        std::shared_ptr<FieldVariable1> fieldVariable1 = (*transferableSolutionData1)[fiberIndex]->variable1[fromVectorIndex].values;
        int componentNo1                               = (*transferableSolutionData1)[fiberIndex]->variable1[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
  }

  // for the second vector of variables (the "algebraics" in case of CellMLAdapter)
  for (int i = 0; i < transferableSolutionData1Front->variable2.size(); i++)
  {
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = slotsConnection.getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b> FieldVariable1;

    std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1Front->variable2[fromVectorIndex].values;
    int componentNo1                               = transferableSolutionData1Front->variable2[fromVectorIndex].componentNo;

    LOG(DEBUG) << "map slot from variable1, index " << fromVectorIndex << " (" << fieldVariable1->name() << "[" << componentNo1 << "])"
      << " to variable" << toVectorNo+1 << ", index " << toVectorIndex;

    if (componentNo1 < 0)
    {
      LOG(DEBUG) << "do not map this slot";
      continue;
    }

    if (toVectorNo == 0)
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable1[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2->variable1[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(7)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData1->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable2 to transferableSolutionData2->variable1

        typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b> FieldVariable1;
        std::shared_ptr<FieldVariable1> fieldVariable1 = (*transferableSolutionData1)[fiberIndex]->variable2[fromVectorIndex].values;
        int componentNo1                               = (*transferableSolutionData1)[fiberIndex]->variable2[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
    else
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable2[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2->variable2[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(8)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData1->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable2 to transferableSolutionData2->variable2

        typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b> FieldVariable1;
        std::shared_ptr<FieldVariable1> fieldVariable1 = (*transferableSolutionData1)[fiberIndex]->variable2[fromVectorIndex].values;
        int componentNo1                               = (*transferableSolutionData1)[fiberIndex]->variable2[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
  }

  // there is no geometry transfer from fibers to anything else
}

/** Transfer between a normal entry and a vector, the first entry of the vector is used
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void SlotConnectorDataTransfer<
  Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  std::vector<std::shared_ptr<
    Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
  >>
>::transfer(const std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
            std::shared_ptr<std::vector<std::shared_ptr<
              Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
            >>> transferableSolutionData2,
            SlotsConnection &slotsConnection,
            int offsetSlotNoData1, int offsetSlotNoData2,
            bool originalTransferDirection)
{
  LOG(DEBUG) << "transfer vector (3)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "offsetSlot: " << offsetSlotNoData1 << ", " << offsetSlotNoData2 << ", slotsConnection: " << slotsConnection.getDebugInformation();
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<std::vector<std::shared_ptr<
          Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
        >>>::getString(transferableSolutionData2);
#endif

  VLOG(1) << "Solution vector mapping (slot_connector_data_transfer_vector.tpp)";

  if (transferableSolutionData2->empty())
  {
    LOG(ERROR) << "Trying to transfer data from a single entry to an empty vector. Types: " << std::endl
      << StringUtility::demangle(typeid(Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>).name()) << std::endl
      << StringUtility::demangle(typeid(std::vector<std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>>>).name());
    return;
  }
/*
  SlotConnectorDataTransfer<SlotConnectorDataType1,SlotConnectorDataType2>::transfer(
    transferableSolutionData1,
    (*transferableSolutionData2)[0],
    slotsConnection
  );*/

  // the following is analogous to the code in slot_connector_data_transfer_fibers_emg.tpp where there are two nested vectors, here it is only one

  assert(transferableSolutionData2->size() > 0);

  std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> transferableSolutionData2Front
    = (*transferableSolutionData2)[0];

  // initialize output connection object
  slotsConnection.initialize(*transferableSolutionData1, *transferableSolutionData2Front, offsetSlotNoData1, offsetSlotNoData2);

  // for the first vector of variables (the "states" in case of CellMLAdapter)
  for (int i = 0; i < transferableSolutionData1->variable1.size(); i++)
  {
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = slotsConnection.getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a> FieldVariable1;

    std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1->variable1[fromVectorIndex].values;
    int componentNo1                               = transferableSolutionData1->variable1[fromVectorIndex].componentNo;

    LOG(DEBUG) << "map slot from variable1, index " << fromVectorIndex << " (" << fieldVariable1->name() << "[" << componentNo1 << "])"
      << " to variable" << toVectorNo+1 << ", index " << toVectorIndex;

    if (componentNo1 < 0)
    {
      LOG(DEBUG) << "do not map this slot";
      continue;
    }

    if (toVectorNo == 0)
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2Front->variable1[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2Front->variable1[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(9)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData2->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable1 to transferableSolutionData2->variable1

        typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
        std::shared_ptr<FieldVariable2> &fieldVariable2 = (*transferableSolutionData2)[fiberIndex]->variable1[fromVectorIndex].values;
        int componentNo2                                = (*transferableSolutionData2)[fiberIndex]->variable1[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
    else
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2Front->variable2[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2Front->variable2[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(10)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData2->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable1 to transferableSolutionData2->variable2

        typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
        std::shared_ptr<FieldVariable2> &fieldVariable2 = (*transferableSolutionData2)[fiberIndex]->variable2[fromVectorIndex].values;
        int componentNo2                                = (*transferableSolutionData2)[fiberIndex]->variable2[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
  }

  // for the second vector of variables (the "algebraics" in case of CellMLAdapter)
  for (int i = 0; i < transferableSolutionData1->variable2.size(); i++)
  {
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = slotsConnection.getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b> FieldVariable1;

    std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1->variable2[fromVectorIndex].values;
    int componentNo1                               = transferableSolutionData1->variable2[fromVectorIndex].componentNo;

    LOG(DEBUG) << "map slot from variable1, index " << fromVectorIndex
      << " to variable" << toVectorNo+1 << ", index " << toVectorIndex;

    if (componentNo1 < 0)
    {
      LOG(DEBUG) << "do not map this slot";
      continue;
    }

    if (toVectorNo == 0)
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2Front->variable1[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2Front->variable1[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(11)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData2->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable2 to transferableSolutionData2->variable1

        typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
        std::shared_ptr<FieldVariable2> &fieldVariable2 = (*transferableSolutionData2)[fiberIndex]->variable1[fromVectorIndex].values;
        int componentNo2                                = (*transferableSolutionData2)[fiberIndex]->variable1[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
    else
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2Front->variable2[toVectorIndex].values;
      int componentNo2                                = transferableSolutionData2Front->variable2[toVectorIndex].componentNo;

      LOG(DEBUG) << "  " << fieldVariable1->name() << "." << fieldVariable1->componentName(componentNo1) << " [" << componentNo1 << "] -> "
        << fieldVariable2->name() << "." << fieldVariable2->componentName(componentNo2) << " [" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible << "(12)";

      // initialize the mapping
      DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo2);

      // loop over vector of fibers
      for (int fiberIndex = 0; fiberIndex < transferableSolutionData2->size(); fiberIndex++)
      {
        // map from transferableSolutionData1->variable2 to transferableSolutionData2->variable2

        typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
        std::shared_ptr<FieldVariable2> &fieldVariable2 = (*transferableSolutionData2)[fiberIndex]->variable2[fromVectorIndex].values;
        int componentNo2                                = (*transferableSolutionData2)[fiberIndex]->variable2[fromVectorIndex].componentNo;

        // map values
        DihuContext::mappingBetweenMeshesManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
      }

      // finalize the mapping
      DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2, componentNo1, componentNo2, avoidCopyIfPossible);
    }
  }

  // transfer geometry field if it was set
  if (transferableSolutionData1->geometryField && !transferableSolutionData2Front->variable1.empty())
  {
    LOG(DEBUG) << "transfer geometry field (1)";

    // get source field variable, this is the same for all fibers
    typedef FieldVariable::FieldVariable<FunctionSpaceType1,3> FieldVariableSource;
    typedef FieldVariable::FieldVariable<FunctionSpaceType2,3> FieldVariableTarget;

    std::shared_ptr<FieldVariableSource> geometryFieldSource = transferableSolutionData1->geometryField;

    // the following prepareMapping would not be necessary for mapping from high to low dimension, as is the case here, mapping from 3D domain to 1D fibers
    std::shared_ptr<FieldVariableTarget> geometryFieldTarget = std::make_shared<FieldVariableTarget>(transferableSolutionData2Front->variable1[0].values->functionSpace()->geometryField());
    DihuContext::mappingBetweenMeshesManager()->template prepareMapping<FieldVariableSource,FieldVariableTarget>(geometryFieldSource, geometryFieldTarget, -1);   // -1 means all components are transferred

    // loop over vector of fibers
    for (int fiberIndex = 0; fiberIndex < transferableSolutionData2->size(); fiberIndex++)
    {
      // map from transferableSolutionData1->geometryField to transferableSolutionData2->variable1[0]->geometryField
      std::shared_ptr<FieldVariableTarget> geometryFieldTarget = std::make_shared<FieldVariableTarget>((*transferableSolutionData2)[fiberIndex]->variable1[0].values->functionSpace()->geometryField());

      // map the whole geometry field (all components, -1), do not avoid copy
      DihuContext::mappingBetweenMeshesManager()->template map<FieldVariableSource,FieldVariableTarget>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
    }

    DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<FieldVariableSource,FieldVariableTarget>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
  }
}
