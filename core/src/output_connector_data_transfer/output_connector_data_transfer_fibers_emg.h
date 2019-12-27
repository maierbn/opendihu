#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "output_connector_data_transfer/output_connection.h"

/** Transfer between the output from cubes partitioned fibers (MultipleInstances<Strang<...) and StaticBidomainSolver
 *
 * template <int nStates, int nIntermediates, typename FunctionSpaceType>
 * struct CellMLOutputConnectorDataType
 * {
 *   Data::FieldVariableComponent<FunctionSpaceType,nStates> stateVariable;          //< one component of the states
 *   Data::FieldVariableComponent<FunctionSpaceType,nIntermediates> intermediateVariable;   //< one component of the intermediates
 * };
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::vector<std::vector<
    Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>   // 1D fibers
  >>,
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>  // 3D field variable
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::vector<std::vector<
                         Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>   // 1D fibers
                       >> &transferableSolutionData1,
                       Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2,
                       OutputConnection &outputConnection);
};

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,  // <3D field variable>
  std::vector<std::vector<
    Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>   // 1D fibers
  >>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                       std::vector<std::vector<
                         Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>   // 1D fibers
                       >> &transferableSolutionData2,
                       OutputConnection &outputConnection);
};

#include "output_connector_data_transfer/output_connector_data_transfer_fibers_emg.tpp"
