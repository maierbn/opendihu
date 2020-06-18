#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <tuple>

#include "output_connector_data_transfer/output_connection.h"

/**
 * Tuple
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2, typename FunctionSpaceType3, int nComponents3a, int nComponents3b>
class SolutionVectorMapping<
  std::tuple<
    std::shared_ptr<OutputConnectorDataType1>,
    std::shared_ptr<OutputConnectorDataType2>
  >,
  Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData3, as efficient as possible
  static void transfer(const std::shared_ptr<std::tuple<
                         std::shared_ptr<OutputConnectorDataType1>,
                         std::shared_ptr<OutputConnectorDataType2>
                       >> transferableSolutionData1,
                       std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType3,nComponents3a,nComponents3b>> transferableSolutionData3,
                       OutputConnection &outputConnection,
                       int offsetSlotNoData1=0, int offsetSlotNoData2=0
                      );
};

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename OutputConnectorDataType2, typename OutputConnectorDataType3>
class SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  std::tuple<
    std::shared_ptr<OutputConnectorDataType2>,
    std::shared_ptr<OutputConnectorDataType3>
  >
>
{
public:

  //! transfer the data from transferableSolutionData1 to transferableSolutionData3, as efficient as possible
  static void transfer(const std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
                       std::shared_ptr<std::tuple<
                         std::shared_ptr<OutputConnectorDataType2>,
                         std::shared_ptr<OutputConnectorDataType3>
                       >> transferableSolutionData2,
                       OutputConnection &outputConnection, int offsetSlotNoData1=0, int offsetSlotNoData2=0);
};

#include "output_connector_data_transfer/output_connector_data_transfer_tuple.tpp"
