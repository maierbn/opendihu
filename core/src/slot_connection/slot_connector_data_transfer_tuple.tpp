#include "slot_connection/slot_connector_data_transfer_tuple.h"

#include "slot_connection/data_helper/slot_connector_data_helper.h"

/** Transfer between a tuple and another field variables with given component number
 */
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2, typename SlotConnectorDataType3>
void SlotConnectorDataTransfer<
  std::tuple<
    std::shared_ptr<SlotConnectorDataType1>,
    std::shared_ptr<SlotConnectorDataType2>
  >,
  SlotConnectorDataType3
>
::transfer(const std::shared_ptr<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
           >> transferableSolutionData1,
           std::shared_ptr<SlotConnectorDataType3> transferableSolutionData3,
           SlotsConnection &slotsConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (1a), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
        >>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<SlotConnectorDataType3>::getString(transferableSolutionData3);
#endif

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection1())
  {
    slotsConnection.subSlotsConnection1() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // do transfer with slots in first tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType3
  >::transfer(std::get<0>(*transferableSolutionData1), transferableSolutionData3, *slotsConnection.subSlotsConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (1b), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<std::tuple<
             std::shared_ptr<SlotConnectorDataType1>,
             std::shared_ptr<SlotConnectorDataType2>
        >>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<SlotConnectorDataType3>::getString(transferableSolutionData3);
#endif

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection2())
  {
    slotsConnection.subSlotsConnection2() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 += SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset does not change
  offsetSlotNoData2 += 0;

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType2,
    SlotConnectorDataType3
  >::transfer(std::get<1>(*transferableSolutionData1), transferableSolutionData3, *slotsConnection.subSlotsConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);
}

/** Transfer between a tuple and another field variables with given component number
 */
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2, typename SlotConnectorDataType3>
void SlotConnectorDataTransfer<
  SlotConnectorDataType1,
  std::tuple<
    std::shared_ptr<SlotConnectorDataType2>,
    std::shared_ptr<SlotConnectorDataType3>
  >
>
::transfer(const std::shared_ptr<SlotConnectorDataType1> transferableSolutionData1,
           std::shared_ptr<std::tuple<
             std::shared_ptr<SlotConnectorDataType2>,
             std::shared_ptr<SlotConnectorDataType3>
           >> transferableSolutionData2,
           SlotsConnection &slotsConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (2a)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<SlotConnectorDataType1>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<std::tuple<
          std::shared_ptr<SlotConnectorDataType2>,
          std::shared_ptr<SlotConnectorDataType3>
        >>::getString(transferableSolutionData2);
#endif

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection1())
  {
    slotsConnection.subSlotsConnection1() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // do transfer with slots in first tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType2
  >::transfer(transferableSolutionData1, std::get<0>(*transferableSolutionData2), *slotsConnection.subSlotsConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (2b)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << SlotConnectorDataHelper<SlotConnectorDataType1>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << SlotConnectorDataHelper<std::tuple<
          std::shared_ptr<SlotConnectorDataType2>,
          std::shared_ptr<SlotConnectorDataType3>
        >>::getString(transferableSolutionData2);
#endif

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection2())
  {
    slotsConnection.subSlotsConnection2() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += SlotConnectorDataHelper<SlotConnectorDataType2>::nSlots(std::get<0>(*transferableSolutionData2));

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType3
  >::transfer(transferableSolutionData1, std::get<1>(*transferableSolutionData2), *slotsConnection.subSlotsConnection2(),
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
           SlotsConnection &slotsConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (3aa)";

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection1())
  {
    slotsConnection.subSlotsConnection1() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  int initialOffsetSlotNoData1 = offsetSlotNoData1;
  int initialOffsetSlotNoData2 = offsetSlotNoData2;

  // do transfer with slots in first tuple part to first tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType3
  >::transfer(std::get<0>(*transferableSolutionData1), std::get<0>(*transferableSolutionData2), *slotsConnection.subSlotsConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3ab)";

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection2())
  {
    slotsConnection.subSlotsConnection2() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += SlotConnectorDataHelper<SlotConnectorDataType3>::nSlots(std::get<0>(*transferableSolutionData2));

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType1,
    SlotConnectorDataType4
  >::transfer(std::get<0>(*transferableSolutionData1), std::get<1>(*transferableSolutionData2), *slotsConnection.subSlotsConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3ba)";
  // reset offsets

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 = initialOffsetSlotNoData1 + SlotConnectorDataHelper<SlotConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset is reset to the start
  offsetSlotNoData2 = initialOffsetSlotNoData2;

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection3())
  {
    slotsConnection.subSlotsConnection3() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType2,
    SlotConnectorDataType3
  >::transfer(std::get<1>(*transferableSolutionData1), std::get<0>(*transferableSolutionData2), *slotsConnection.subSlotsConnection3(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3bb)";

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += SlotConnectorDataHelper<SlotConnectorDataType3>::nSlots(std::get<0>(*transferableSolutionData2));

  // copy slotsConnection for second tuple
  if (!slotsConnection.subSlotsConnection4())
  {
    slotsConnection.subSlotsConnection4() = std::make_shared<SlotsConnection>(slotsConnection);
  }

  // do transfer with slots in second tuple part
  SlotConnectorDataTransfer<
    SlotConnectorDataType2,
    SlotConnectorDataType4
  >::transfer(std::get<1>(*transferableSolutionData1), std::get<1>(*transferableSolutionData2), *slotsConnection.subSlotsConnection4(),
              offsetSlotNoData1, offsetSlotNoData2);
}
