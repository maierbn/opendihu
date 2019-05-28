#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

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
template<typename TransferableSolutionDataType1, typename TransferableSolutionDataType2>
class SolutionVectorMapping
{
};

/** Transfer between two field variables with given component number, both field variables have a component no. != 1
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2>
class SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,   // <fieldVariableType,componentNo,prefactor>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>,int,double> &transferableSolutionData1,
                       const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double> &transferableSolutionData2);
};

/** Transfer between two field variables with given component number,
 *  the second field variable has only 1 component
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2>
class SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,   // <fieldVariableType,componentNo,prefactor>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>, int, double>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>,int,double> &transferableSolutionData1,
                       const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>,int,double> &transferableSolutionData2);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.tpp"

// include all other solution_vector_mapping partial specializations
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_finite_element_method_cellml.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_multidomain.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_vector.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_bidomain.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_cellml.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_fibers_emg.h"
