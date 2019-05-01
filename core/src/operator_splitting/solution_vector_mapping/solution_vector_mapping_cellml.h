#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

/** Transfer between a pair of states and intermediates field variables from CellML, with given component numbers each, to a normal field variable with component no, both field variables have a component no. != 1
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2>
class SolutionVectorMapping<
  std::pair<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
  >,
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::pair<
                         std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
                         std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
                       > &transferableSolutionData1,
                       const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double> &transferableSolutionData2);
};

/** Reverse transfer to the one before
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,
  std::pair<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
  >
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>,int,double> &transferableSolutionData1,
                      const std::pair<
                        std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
                        std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
                      > &transferableSolutionData2);
};

/** Transfer between a pair of states and intermediates field variables from CellML, with given component numbers each, to a normal field variable with component no,
 *  the second field variable has only 1 component
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2>
class SolutionVectorMapping<
  std::pair<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
  >,
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>, int, double>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::pair<
                         std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
                         std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
                       > &transferableSolutionData1,
                       const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>,int,double> &transferableSolutionData2);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_cellml.tpp"
