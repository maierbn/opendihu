#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

/** Transfer between the output from cubes partitioned fibers (MultipleInstances<Strang<...) and StaticBidomainSolver
 *
 * template <int nStates, int nIntermediates, typename FunctionSpaceType>
 * struct CellMLOutputConnectorDataType
 * {
 *   Data::ScaledFieldVariableComponent<FunctionSpaceType,nStates> stateVariable;          //< one component of the states
 *   Data::ScaledFieldVariableComponent<FunctionSpaceType,nIntermediates> intermediateVariable;   //< one component of the intermediates
 * };
 */
template<typename BasisFunctionType, int nComponents1a, int nComponents1b, typename FieldVariableType2>
class SolutionVectorMapping<
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents1a,nComponents1b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>,
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::vector<std::vector<
                         CellMLOutputConnectorDataType<
                           nComponents1a,nComponents1b,
                           FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
                         >
                       >> &transferableSolutionData1,
                       std::shared_ptr<FieldVariableType2> transferableSolutionData2,
                       const std::string transferSlotName);
};

template<typename BasisFunctionType, typename FieldVariableType1, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents2a,nComponents2b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(std::shared_ptr<FieldVariableType1> transferableSolutionData1,
                       const std::vector<std::vector<
                         CellMLOutputConnectorDataType<
                           nComponents2a,nComponents2b,
                           FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
                         >
                       >> &transferableSolutionData2,
                       const std::string transferSlotName);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_fibers_emg.tpp"
