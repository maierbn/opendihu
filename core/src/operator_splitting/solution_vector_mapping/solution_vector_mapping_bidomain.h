#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

/** Transfer between the output from cubes partitioned fibers (MultipleInstances<Strang<...) and StaticBidomainSolver
 */
template<typename BasisFunctionType, int nComponents1, typename FieldVariableType2>
class SolutionVectorMapping<
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1>>,
      int, double>
    >>,   // vector<vector<fieldVariableType,componentNo,prefactor>>
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::vector<std::vector<
                       std::tuple<
                         std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1>>,
                         int, double>
                       >> &transferableSolutionData1,
                       std::shared_ptr<FieldVariableType2> transferableSolutionData2,
                       const std::string transferSlotName);
};

template<typename BasisFunctionType, typename FieldVariableType1, int nComponents2>
class SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents2>>,
      int, double>
    >>   // vector<vector<fieldVariableType,componentNo,prefactor>>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(std::shared_ptr<FieldVariableType1> transferableSolutionData1,
                       std::vector<std::vector<
                       std::tuple<
                         std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents2>>,
                         int, double>
                       >> transferableSolutionData2,
                       const std::string transferSlotName);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_bidomain.tpp"

/*
/store/software/opendihu/core/src/operator_splitting/coupling_or_godunov.tpp:53:15: error: ‘transfer’ is not a member of ‘
SolutionVectorMapping<
  std::vector<
    std::vector<
      std::tuple<
        std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >,
        int,
        double
      >,
      std::allocator<
        std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >,
          int, double>
      >
    >,
    std::allocator<
      std::vector<
        std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >, int, double
        >,
        std::allocator<
          std::tuple<
            std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >, int, double
          >
        >
      >
    >
  >,

  std::pair<
    std::vector<_p_Vec*>,
    std::vector<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1> >, 1> >, std::allocator<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1> >, 1> > >
    >
  >
>
’
 */
