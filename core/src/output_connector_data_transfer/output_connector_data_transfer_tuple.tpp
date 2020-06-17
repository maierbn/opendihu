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

  // do transfer with slots in first tuple part
  SolutionVectorMapping<
    OutputConnectorDataType1,
    Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
  >::transfer(std::get<0>(*transferableSolutionData1), transferableSolutionData3, outputConnection,
              offsetSlotNoData1, offsetSlotNoData2);

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection())
  {
    outputConnection.subOutputConnection() = std::make_shared<OutputConnection>(outputConnection);
  }

  // first offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData1 += OutputConnectorDataHelper<OutputConnectorDataType1>::nSlots(std::get<0>(*transferableSolutionData1));

  // second offset does not change
  offsetSlotNoData2 += 0;

  LOG(DEBUG) << "transfer tuple (1b), offsetSlotNoData1: " << offsetSlotNoData1 << ", offsetSlotNoData2: " << offsetSlotNoData2;

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    OutputConnectorDataType2,
    Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
  >::transfer(std::get<1>(*transferableSolutionData1), transferableSolutionData3, *outputConnection.subOutputConnection(),
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

  // do transfer with slots in first tuple part
  SolutionVectorMapping<
    Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
    OutputConnectorDataType2
  >::transfer(transferableSolutionData1, std::get<0>(*transferableSolutionData2), outputConnection,
              offsetSlotNoData1, offsetSlotNoData2);

  // copy outputConnection for second tuple
  if (!outputConnection.subOutputConnection())
  {
    outputConnection.subOutputConnection() = std::make_shared<OutputConnection>(outputConnection);
  }

  // first offset does not change
  offsetSlotNoData1 += 0;

  // second offset changes by number of slots that were present in the first tuple entry
  offsetSlotNoData2 += OutputConnectorDataHelper<OutputConnectorDataType2>::nSlots(std::get<0>(*transferableSolutionData2));

  LOG(DEBUG) << "transfer tuple (2b)";

  // do transfer with slots in second tuple part
  SolutionVectorMapping<
    Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
    OutputConnectorDataType3
  >::transfer(transferableSolutionData1, std::get<1>(*transferableSolutionData2), *outputConnection.subOutputConnection(),
              offsetSlotNoData1, offsetSlotNoData2);
}
