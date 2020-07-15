#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <tuple>

#include "slot_connection/slots_connection.h"

/**
 * Tuple
 */
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2, typename FunctionSpaceType3, int nComponents3a, int nComponents3b>
class SlotConnectorDataTransfer<
  std::tuple<
    std::shared_ptr<SlotConnectorDataType1>,
    std::shared_ptr<SlotConnectorDataType2>
  >,
  Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData3, as efficient as possible
  static void transfer(const std::shared_ptr<std::tuple<
                         std::shared_ptr<SlotConnectorDataType1>,
                         std::shared_ptr<SlotConnectorDataType2>
                       >> transferableSolutionData1,
                       std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>> transferableSolutionData3,
                       SlotsConnection &slotsConnection,
                       int offsetSlotNoData1=0, int offsetSlotNoData2=0
                      );
};

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename SlotConnectorDataType2, typename SlotConnectorDataType3>
class SlotConnectorDataTransfer<
  Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  std::tuple<
    std::shared_ptr<SlotConnectorDataType2>,
    std::shared_ptr<SlotConnectorDataType3>
  >
>
{
public:

  //! transfer the data from transferableSolutionData1 to transferableSolutionData3, as efficient as possible
  static void transfer(const std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
                       std::shared_ptr<std::tuple<
                         std::shared_ptr<SlotConnectorDataType2>,
                         std::shared_ptr<SlotConnectorDataType3>
                       >> transferableSolutionData2,
                       SlotsConnection &slotsConnection, int offsetSlotNoData1=0, int offsetSlotNoData2=0);
};


template<typename SlotConnectorDataType1, typename SlotConnectorDataType2, typename SlotConnectorDataType3, typename SlotConnectorDataType4>
class SlotConnectorDataTransfer<
  std::tuple<
    std::shared_ptr<SlotConnectorDataType1>,
    std::shared_ptr<SlotConnectorDataType2>
  >,
  std::tuple<
    std::shared_ptr<SlotConnectorDataType3>,
    std::shared_ptr<SlotConnectorDataType4>
  >
>
{
public:

  //! transfer the data from transferableSolutionData1 to transferableSolutionData3, as efficient as possible
  static void transfer(const std::shared_ptr<std::tuple<
                         std::shared_ptr<SlotConnectorDataType1>,
                         std::shared_ptr<SlotConnectorDataType2>
                       >> transferableSolutionData1,
                       std::shared_ptr<std::tuple<
                         std::shared_ptr<SlotConnectorDataType3>,
                         std::shared_ptr<SlotConnectorDataType4>
                       >> transferableSolutionData2,
                       SlotsConnection &slotsConnection, int offsetSlotNoData1=0, int offsetSlotNoData2=0);
};

#include "slot_connection/slot_connector_data_transfer_tuple.tpp"
