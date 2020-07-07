#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "output_connector_data_transfer/output_connection.h"

/** Transfer between the output from cubes partitioned fibers (MultipleInstances<Strang<MultipleInstances<...) and StaticBidomainSolver
 *
 * template <int nStates, int nAlgebraics, typename FunctionSpaceType>
 * struct CellMLOutputConnectorDataType
 * {
 *   Data::FieldVariableComponent<FunctionSpaceType,nStates> stateVariable;          //< one component of the states
 *   Data::FieldVariableComponent<FunctionSpaceType,nAlgebraics> algebraicVariable;   //< one component of the algebraics
 * };
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename OutputConnectorDataType1c, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::vector<std::shared_ptr<
    std::tuple<std::shared_ptr<
      std::vector<std::shared_ptr<
        Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>   // 1D fibers
      >>>,
      OutputConnectorDataType1c
    >
  >>,
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>  // 3D field variable
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<std::vector<std::shared_ptr<
                         std::tuple<std::shared_ptr<
                           std::vector<std::shared_ptr<
                             Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>   // 1D fibers
                           >>>,
                           OutputConnectorDataType1c
                         >
                       >>> transferableSolutionData1,
                       std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> transferableSolutionData2,
                       OutputConnection &outputConnection,
                       int offsetSlotNoData1=0, int offsetSlotNoData2=0);
};

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b, typename OutputConnectorDataType2c>
class SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,  // <3D field variable>
  std::vector<std::shared_ptr<
    std::tuple<std::shared_ptr<
      std::vector<std::shared_ptr<
        Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>   // 1D fibers
      >>>,
      OutputConnectorDataType2c
    >
  >>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
                       std::shared_ptr<std::vector<std::shared_ptr<
                          std::tuple<std::shared_ptr<
                            std::vector<std::shared_ptr<
                              Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>   // 1D fibers
                            >>>,
                            OutputConnectorDataType2c
                          >
                        >>> transferableSolutionData2,
                       OutputConnection &outputConnection,
                       int offsetSlotNoData1=0, int offsetSlotNoData2=0);
};

#include "output_connector_data_transfer/output_connector_data_transfer_fibers_emg.tpp"
