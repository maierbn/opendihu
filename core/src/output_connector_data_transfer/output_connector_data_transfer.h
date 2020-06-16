#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/00_cellml_adapter_base.h"
#include "output_connector_data_transfer/output_connection.h"

/**
 * The Data classes contain each a vector that stores the solution. Often, the values need to be accessed to
 * continue computation using them, e.g. in operator splittings. Then possibly not all the values need to be accessed
 * but only some of them while the others are only needed for the own computation. An example is a cellml model that
 * contains lots of internal states and only a single state variable is needed in the diffusion equation.
 *
 * Because the data may be stored in different forms, e.g. as raw Vec, as FieldVariable with several components where only one component is the solution variable
 * or as nested vector, the transfer function that transfers the data from one object to another is implemented using partial specializations of this class.
 *
 * The data is also scaled with a prefactor, that can be given. This is needed for the results of the cellml class.
 */
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
class SolutionVectorMapping
{
};

/** Transfer between a pair of states and algebraics field variables from CellML, with given component numbers each, to a normal field variable with component no, both field variables have a component no. != 1
 *
 *   template<typename FunctionSpaceType, int nComponents1, int nComponents2>
 *   struct OutputConnectorData
 *   {
 *     std::vector<ComponentOfFieldVariable<FunctionSpaceType,nComponents1>> variable1;    //< vector of indications of components of field variables, the field variables have all the same number of components
 *     std::vector<ComponentOfFieldVariable<FunctionSpaceType,nComponents2>> variable2;    //< second vector with different number of components for the field variables
 *   }
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,      //< algebraics and states from cellmlAdapter
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
                       std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> transferableSolutionData2,
                       OutputConnection &outputConnection,
                       int offsetSlotNoData1=0, int offsetSlotNoData2=0
                      );
};

#include "output_connector_data_transfer/output_connector_data_transfer.tpp"
#include "output_connector_data_transfer/output_connector_data_transfer_fibers_emg.h"
#include "output_connector_data_transfer/output_connector_data_transfer_vector.h"
#include "output_connector_data_transfer/output_connector_data_transfer_tuple.h"
