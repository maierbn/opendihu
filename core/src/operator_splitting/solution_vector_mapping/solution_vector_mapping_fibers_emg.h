#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

/** Transfer between the output from cubes partitioned fibers (MultipleInstances<Strang<...) and StaticBidomainSolver
 */
template<typename BasisFunctionType, int nComponents1a, int nComponents1b, typename FieldVariableType2>
class SolutionVectorMapping<
  std::vector<std::vector<
    std::pair<
      std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents1a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
      std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents1b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
    >
  >>,
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(const std::vector<std::vector<
                         std::pair<
                           std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents1a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
                           std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents1b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
                         >
                       >> &transferableSolutionData1,
                       std::shared_ptr<FieldVariableType2> transferableSolutionData2);
};

template<typename BasisFunctionType, typename FieldVariableType1, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    std::pair<
      std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents2a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
      std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents2b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
    >
  >>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible
  static void transfer(std::shared_ptr<FieldVariableType1> transferableSolutionData1,
                       const std::vector<std::vector<
                         std::pair<
                           std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents2a>>, int, double>,   // <fieldVariableTypeStates,componentNoStates,prefactor>
                           std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>,nComponents2b>>, int>    // <fieldIariableIntermediates,componentNoIntermediates
                         >
                       >> &transferableSolutionData2);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_fibers_emg.tpp"
