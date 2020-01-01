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

/** Transfer between a pair of states and intermediates field variables from CellML, with given component numbers each, to a normal field variable with component no, both field variables have a component no. != 1
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
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,      //< intermediates and states from cellmlAdapter
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                       Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2,
                       OutputConnection &outputConnection);
};

/** Reverse transfer to the one before
 */
/*
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1>,
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
                       const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2,
                       const std::string transferSlotName);
};
*/

/** Transfer between a pair of states and intermediates field variables from CellML, with given component numbers each, to a normal field variable with component no,
 *  the second field variable has only 1 component
 */
/*template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2>
class SolutionVectorMapping<
  OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  Data::FieldVariableComponent<FunctionSpaceType2,1>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                       const Data::FieldVariableComponent<FunctionSpaceType2,1> &transferableSolutionData2,
                       const std::string transferSlotName);
};*/

#include "output_connector_data_transfer/output_connector_data_transfer_cellml.tpp"
