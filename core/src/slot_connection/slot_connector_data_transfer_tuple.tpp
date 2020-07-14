#include "slot_connection/slot_connector_data_transfer_tuple.h"

#include "slot_connection/slot_connector_data_helper.h"

/** Transfer between a tuple and another field variables with given component number
 */
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2, typename FunctionSpaceType3, int nComponents3a, int nComponents3b>
void SlotConnectorDataTransfer<
  std::tuple<
    std::shared_ptr<SlotConnectorDataType1>,
    std::shared_ptr<SlotConnectorDataType2>
  >,
  Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
>
::transfer(const std::shared_ptr<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
           >> transferableSolutionData1,
           std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>> transferableSolutionData3,
           SlotConnection &outputConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (1a), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
        >>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>>::getString(transferableSolutionData3);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection1())
  {
    outputConnection.subSlotConnection1() = std::make_shared<SlotConnection>(outputConnection);
  }

  // do transfer with slots in first tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
  >::transfer(std::get<0>(*transferableSolutionData1), transferableSolutionData3, *outputConnection.subSlotConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (1b), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
        >>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>>::getString(transferableSolutionData3);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection2())
  {
    outputConnection.subSlotConnection2() = std::make_shared<SlotConnection>(outputConnection);
  }

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 += SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset does not change
  offsetSlotNoData2 += 0;

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType2,
    Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
  >::transfer(std::get<1>(*transferableSolutionData1), transferableSolutionData3, *outputConnection.subSlotConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);
}

/** Transfer between a tuple and another field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename SlotConnectorDataType2, typename SlotConnectorDataType3>
void SlotConnectorDataTransfer<
  Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  std::tuple<
    std::shared_ptr<SlotConnectorDataType2>,
    std::shared_ptr<SlotConnectorDataType3>
  >
>
::transfer(const std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
           std::shared_ptr<std::tuple<
             std::shared_ptr<SlotConnectorDataType2>,
             std::shared_ptr<SlotConnectorDataType3>
           >> transferableSolutionData2,
           SlotConnection &outputConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (2a)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<std::tuple<
          std::shared_ptr<SlotConnectorDataType2>,
          std::shared_ptr<SlotConnectorDataType3>
        >>::getString(transferableSolutionData2);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection1())
  {
    outputConnection.subSlotConnection1() = std::make_shared<SlotConnection>(outputConnection);
  }

  // do transfer with slots in first tuple part
  SlotConnectorDataTransfer<
    Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
    SlotConnectorDataType2
  >::transfer(transferableSolutionData1, std::get<0>(*transferableSolutionData2), *outputConnection.subSlotConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (2b)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<std::tuple<
          std::shared_ptr<SlotConnectorDataType2>,
          std::shared_ptr<SlotConnectorDataType3>
        >>::getString(transferableSolutionData2);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection2())
  {
    outputConnection.subSlotConnection2() = std::make_shared<SlotConnection>(outputConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<0>(*transferableSolutionData2));

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
    SlotConnectorDataType3
  >::transfer(transferableSolutionData1, std::get<1>(*transferableSolutionData2), *outputConnection.subSlotConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);
}

template<typename SlotConnectorDataType1, typename SlotConnectorDataType2, typename SlotConnectorDataType3, typename SlotConnectorDataType4>
void SlotConnectorDataTransfer<
  std::tuple<
    std::shared_ptr<SlotConnectorDataType1>,
    std::shared_ptr<SlotConnectorDataType2>
  >,
  std::tuple<
    std::shared_ptr<SlotConnectorDataType3>,
    std::shared_ptr<SlotConnectorDataType4>
  >
>
::transfer(const std::shared_ptr<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
           >> transferableSolutionData1,
             std::shared_ptr<std::tuple<
             std::shared_ptr<SlotConnectorDataType3>,
             std::shared_ptr<SlotConnectorDataType4>
           >> transferableSolutionData2,
           SlotConnection &outputConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (3aa)";

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection1())
  {
    outputConnection.subSlotConnection1() = std::make_shared<SlotConnection>(outputConnection);
  }

  int initialOffsetSlotNoData1 = offsetSlotNoData1;
  int initialOffsetSlotNoData2 = offsetSlotNoData2;

  // do transfer with slots in first tuple part to first tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType3
  >::transfer(std::get<0>(*transferableSolutionData1), std::get<0>(*transferableSolutionData2), *outputConnection.subSlotConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3ab)";

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection2())
  {
    outputConnection.subSlotConnection2() = std::make_shared<SlotConnection>(outputConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += SlotConnectorDataHelper<SlotConnectorDataType3>::nSlots(std::get<0>(*transferableSolutionData2));

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType4
  >::transfer(std::get<0>(*transferableSolutionData1), std::get<1>(*transferableSolutionData2), *outputConnection.subSlotConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3ba)";
  // reset offsets

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 = initialOffsetSlotNoData1 + SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset is reset to the start
  offsetSlotNoData2 = initialOffsetSlotNoData2;

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection3())
  {
    outputConnection.subSlotConnection3() = std::make_shared<SlotConnection>(outputConnection);
  }

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType2,
    SlotConnectorDataType3
  >::transfer(std::get<1>(*transferableSolutionData1), std::get<0>(*transferableSolutionData2), *outputConnection.subSlotConnection3(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3bb)";

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += SlotConnectorDataHelper<SlotConnectorDataType3>::nSlots(std::get<0>(*transferableSolutionData2));

  // copy outputConnection for second tuple
  if (!outputConnection.subSlotConnection4())
  {
    outputConnection.subSlotConnection4() = std::make_shared<SlotConnection>(outputConnection);
  }

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType2,
    SlotConnectorDataType4
  >::transfer(std::get<1>(*transferableSolutionData1), std::get<1>(*transferableSolutionData2), *outputConnection.subSlotConnection4(),
              offsetSlotNoData1, offsetSlotNoData2);
}
