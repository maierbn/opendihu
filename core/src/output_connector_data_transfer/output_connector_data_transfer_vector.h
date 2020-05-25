#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "output_connector_data_transfer/output_connection.h"

/**
 * The Data classes contain each a vector that stores the solution. Often, the values need to be accessed to
 * continue computation using them, e.g. in operator splittings. Then possibly not all the values need to be accessed
 * but only some of them while the others are only needed for the own computation. An example is a cellml model that
 * contains lots of internal states and only a single state variable is needed in the diffusion equation.
 *
 * Because the data may be stored in different forms, e.g. as raw Vec, as FieldVariable with several components where only one component is the solution variable
 * or as nested vector, the transfer function that transfers the data from one object to another is implemented using partial specializations of this class.
 */

/** Transfer between two vectors of any type, this calls the implementations without vectors.
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
class SolutionVectorMapping<
  std::vector<std::shared_ptr<OutputConnectorDataType1>>,
  std::vector<std::shared_ptr<OutputConnectorDataType2>>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType1>>> transferableSolutionData1,
                       std::shared_ptr<std::vector<std::shared_ptr<OutputConnectorDataType2>>> transferableSolutionData2,
                       OutputConnection &outputConnection);
};

/** Transfer between a normal entry and a vector, the first entry of the vector is used
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::vector<std::shared_ptr<
    Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>
  >>,
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<std::vector<std::shared_ptr<
                         Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>
                       >>> transferableSolutionData1,
                       std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> transferableSolutionData2,
                       OutputConnection &outputConnection);
};

/** Transfer between a normal entry and a vector, the first entry of the vector is used
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  std::vector<std::shared_ptr<
    Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
  >>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
                       std::shared_ptr<std::vector<std::shared_ptr<
                         Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
                       >>> transferableSolutionData2,
                       OutputConnection &outputConnection);
};

#include "output_connector_data_transfer/output_connector_data_transfer_vector.tpp"
