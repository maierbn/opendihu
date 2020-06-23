#include "output_connector_data_transfer/output_connector_data_transfer_tuple.h"

#include "output_connector_data_transfer/output_connector_data_helper.h"

/** Transfer between a tuple and another field variables with given component number
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2, typename FunctionSpaceType3, int nComponents3a, int nComponents3b>
void SolutionVectorMapping<
  std::tuple<
    std::shared_ptr<OutputConnectorDataType1>,
    std::shared_ptr<OutputConnectorDataType2>
  >,
  Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
>
::transfer(const std::shared_ptr<std::tuple<
             std::shared_ptr<OutputConnectorDataType1>,
             std::shared_ptr<OutputConnectorDataType2>
           >> transferableSolutionData1,
           std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>> transferableSolutionData3,
           OutputConnection &outputConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (1a), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << OutputConnectorDataHelper<std::tuple<
             std::shared_ptr<OutputConnectorDataType1>,
             std::shared_ptr<OutputConnectorDataType2>
        >>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << OutputConnectorDataHelper<Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>>::getString(transferableSolutionData3);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection1())
  {
    outputConnection.subOutputConnection1() = std::make_shared<OutputConnection>(outputConnection);
  }

  // do transfer with slots in first tuple part
  SolutionVectorMapping<
    OutputConnectorDataType1,
    Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
  >::transfer(std::get<0>(*transferableSolutionData1), transferableSolutionData3, *outputConnection.subOutputConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (1b), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << OutputConnectorDataHelper<std::tuple<
             std::shared_ptr<OutputConnectorDataType1>,
             std::shared_ptr<OutputConnectorDataType2>
        >>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << OutputConnectorDataHelper<Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>>::getString(transferableSolutionData3);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection2())
  {
    outputConnection.subOutputConnection2() = std::make_shared<OutputConnection>(outputConnection);
  }

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 += OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset does not change
  offsetSlotNoData2 += 0;

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    OutputConnectorDataType2,
    Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
  >::transfer(std::get<1>(*transferableSolutionData1), transferableSolutionData3, *outputConnection.subOutputConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);
}

/** Transfer between a tuple and another field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename OutputConnectorDataType2, typename OutputConnectorDataType3>
void SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  std::tuple<
    std::shared_ptr<OutputConnectorDataType2>,
    std::shared_ptr<OutputConnectorDataType3>
  >
>
::transfer(const std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
           std::shared_ptr<std::tuple<
             std::shared_ptr<OutputConnectorDataType2>,
             std::shared_ptr<OutputConnectorDataType3>
           >> transferableSolutionData2,
           OutputConnection &outputConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (2a)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << OutputConnectorDataHelper<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << OutputConnectorDataHelper<std::tuple<
          std::shared_ptr<OutputConnectorDataType2>,
          std::shared_ptr<OutputConnectorDataType3>
        >>::getString(transferableSolutionData2);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection1())
  {
    outputConnection.subOutputConnection1() = std::make_shared<OutputConnection>(outputConnection);
  }

  // do transfer with slots in first tuple part
  SolutionVectorMapping<
    Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
    OutputConnectorDataType2
  >::transfer(transferableSolutionData1, std::get<0>(*transferableSolutionData2), *outputConnection.subOutputConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (2b)";
#ifdef SOLUTION_VECTOR_MAPPING_DEBUGGING_OUTPUT
  LOG(DEBUG) << "transferableSolutionData1: "
    << OutputConnectorDataHelper<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>>::getString(transferableSolutionData1);
  LOG(DEBUG) << "transferableSolutionData2: "
    << OutputConnectorDataHelper<std::tuple<
          std::shared_ptr<OutputConnectorDataType2>,
          std::shared_ptr<OutputConnectorDataType3>
        >>::getString(transferableSolutionData2);
#endif

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection2())
  {
    outputConnection.subOutputConnection2() = std::make_shared<OutputConnection>(outputConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<0>(*transferableSolutionData2));

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
    OutputConnectorDataType3
  >::transfer(transferableSolutionData1, std::get<1>(*transferableSolutionData2), *outputConnection.subOutputConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);
}

template<typename OutputConnectorDataType1, typename OutputConnectorDataType2, typename OutputConnectorDataType3, typename OutputConnectorDataType4>
void SolutionVectorMapping<
  std::tuple<
    std::shared_ptr<OutputConnectorDataType1>,
    std::shared_ptr<OutputConnectorDataType2>
  >,
  std::tuple<
    std::shared_ptr<OutputConnectorDataType3>,
    std::shared_ptr<OutputConnectorDataType4>
  >
>
::transfer(const std::shared_ptr<std::tuple<
             std::shared_ptr<OutputConnectorDataType1>,
             std::shared_ptr<OutputConnectorDataType2>
           >> transferableSolutionData1,
             std::shared_ptr<std::tuple<
             std::shared_ptr<OutputConnectorDataType3>,
             std::shared_ptr<OutputConnectorDataType4>
           >> transferableSolutionData2,
           OutputConnection &outputConnection, int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "transfer tuple (3aa)";

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection1())
  {
    outputConnection.subOutputConnection1() = std::make_shared<OutputConnection>(outputConnection);
  }

  int initialOffsetSlotNoData1 = offsetSlotNoData1;
  int initialOffsetSlotNoData2 = offsetSlotNoData2;

  // do transfer with slots in first tuple part to first tuple part
  SolutionVectorMapping<
    OutputConnectorDataType1,
    OutputConnectorDataType3
  >::transfer(std::get<0>(*transferableSolutionData1), std::get<0>(*transferableSolutionData2), *outputConnection.subOutputConnection1(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3ab)";

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection2())
  {
    outputConnection.subOutputConnection2() = std::make_shared<OutputConnection>(outputConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += OutputConnectorDataHelper<OutputConnectorDataType3>::nSlots(std::get<0>(*transferableSolutionData2));

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    OutputConnectorDataType1,
    OutputConnectorDataType4
  >::transfer(std::get<0>(*transferableSolutionData1), std::get<1>(*transferableSolutionData2), *outputConnection.subOutputConnection2(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3ba)";
  // reset offsets

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 = initialOffsetSlotNoData1 + OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset is reset to the start
  offsetSlotNoData2 = initialOffsetSlotNoData2;

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection3())
  {
    outputConnection.subOutputConnection3() = std::make_shared<OutputConnection>(outputConnection);
  }

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    OutputConnectorDataType2,
    OutputConnectorDataType3
  >::transfer(std::get<1>(*transferableSolutionData1), std::get<0>(*transferableSolutionData2), *outputConnection.subOutputConnection3(),
              offsetSlotNoData1, offsetSlotNoData2);

  LOG(DEBUG) << "transfer tuple (3bb)";

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += OutputConnectorDataHelper<OutputConnectorDataType3>::nSlots(std::get<0>(*transferableSolutionData2));

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection4())
  {
    outputConnection.subOutputConnection4() = std::make_shared<OutputConnection>(outputConnection);
  }

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    OutputConnectorDataType2,
    OutputConnectorDataType4
  >::transfer(std::get<1>(*transferableSolutionData1), std::get<1>(*transferableSolutionData2), *outputConnection.subOutputConnection4(),
              offsetSlotNoData1, offsetSlotNoData2);
}
