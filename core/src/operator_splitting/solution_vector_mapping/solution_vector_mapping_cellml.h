#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/00_cellml_adapter_base.h"

/** Transfer between a pair of states and intermediates field variables from CellML, with given component numbers each, to a normal field variable with component no, both field variables have a component no. != 1
 *
 * template <int nStates, int nIntermediates, typename FunctionSpaceType>
 * struct CellMLOutputConnectorDataType
 * {
 *   Data::ScaledFieldVariableComponent<FunctionSpaceType,nStates> stateVariable;          //< one component of the states
 *   Data::ScaledFieldVariableComponent<FunctionSpaceType,nIntermediates> intermediateVariable;   //< one component of the intermediates
 * };
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2>
class SolutionVectorMapping<
  CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1>,      //< intermediates and states from cellmlAdapter
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1> &transferableSolutionData1,
                       const Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2> &transferableSolutionData2,
                       const std::string transferSlotName);
};

/** Reverse transfer to the one before
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>,
  CellMLOutputConnectorDataType<nComponents2a,nComponents2b,FunctionSpaceType2>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
                       const CellMLOutputConnectorDataType<nComponents2a,nComponents2b,FunctionSpaceType2> &transferableSolutionData2,
                       const std::string transferSlotName);
};

/** Transfer between a pair of states and intermediates field variables from CellML, with given component numbers each, to a normal field variable with component no,
 *  the second field variable has only 1 component
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2>
class SolutionVectorMapping<
  CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1>,
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,1>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1> &transferableSolutionData1,
                       const Data::ScaledFieldVariableComponent<FunctionSpaceType2,1> &transferableSolutionData2,
                       const std::string transferSlotName);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_cellml.tpp"
