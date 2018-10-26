#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

/** Transfer between one scalar field variable (e.g. from finite element method) and one field variables with given component number and prefactor (e.g. from cellml)
 */
template<typename FunctionSpaceType1, typename FunctionSpaceType2, int nComponents2>
class SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,   // <fieldVariableType,componentNo>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double>
> :
  public SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,int,double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double>
  >
{
public:
  ///! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
                       const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double> &transferableSolutionData2);

};

/** Transfer between one scalar field variable (e.g. from finite element method) and one field variables with given component number and prefactor (e.g. from cellml)
 *  The other way
 */
template<typename FunctionSpaceType1, typename FunctionSpaceType2, int nComponents2>
class SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double>,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>    // <fieldVariableType,componentNo>
> :
  public SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,int,double>
  >
{
public:
  ///! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>,int,double> &transferableSolutionData2,
                       const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1);

};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_finite_element_method_cellml.tpp"
