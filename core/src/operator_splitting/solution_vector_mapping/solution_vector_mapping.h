#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/time_stepping/time_stepping.h"

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

/** Transfer between two field variables with given component number, both field variables have a component no. != 1
 *
 *  template<typename FunctionSpaceType, int nComponents>
 *  struct ScaledFieldVariableComponent
 *  {
 *    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> values;    //< a field variable containing the payload data that is to be exchangend to another solver
 *    int componentNo;                              //< the component of values that is relevant, only this component out of the potentially multi-component field variable in values will be transferred.
 *    double scalingFactor;    //< a scaling factor, the values will be multiplied by the factor before the transfer. Disabled if set to 1.0.
 *  };
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2>
class SolutionVectorMapping<
  Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>,   // <fieldVariableType,componentNo,prefactor>
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
                       const Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2> &transferableSolutionData2,
                       const std::string transferSlotName
                      );
};

/** Transfer between two field variables with given component number,
 *  the second field variable has only 1 component
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2>
class SolutionVectorMapping<
  Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>,   // <fieldVariableType,componentNo,prefactor>
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,1>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
                       const Data::ScaledFieldVariableComponent<FunctionSpaceType2,1> &transferableSolutionData2,
                       const std::string transferSlotName);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.tpp"

// include all other solution_vector_mapping partial specializations
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_finite_element_method_cellml.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_multidomain.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_vector.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_bidomain.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_cellml.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_fibers_emg.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_linear_elasticity.h"
