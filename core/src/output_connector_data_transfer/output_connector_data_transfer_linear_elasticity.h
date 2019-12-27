#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "cellml/00_cellml_adapter_base.h"
#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_linear_elasticity_solver.h"

/** Transfer between the output from cubes partitioned fibers (MultipleInstances<Strang<...) and an elasticity solver, e.g. QuasiStaticNonlinearElasticitySolverFebio
 *
 * template <int nStates, int nIntermediates, typename FunctionSpaceType>
 * struct CellMLOutputConnectorDataType
 * {
 *   Data::FieldVariableComponent<FunctionSpaceType,nStates> stateVariable;          //< one component of the states
 *   Data::FieldVariableComponent<FunctionSpaceType,nIntermediates> intermediateVariable;   //< one component of the intermediates
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
  TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType2>
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
                       TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType2> transferableSolutionData2,
                       const std::string transferSlotName);
};

template<typename FieldVariableType1, typename BasisFunctionType, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType1>,  // 3D data
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
  static void transfer(TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType1> transferableSolutionData1,
                       const std::vector<std::vector<
                         CellMLOutputConnectorDataType<
                           nComponents2a,nComponents2b,
                           FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
                         >
                       >> &transferableSolutionData2,
                       const std::string transferSlotName);
};

#include "output_connector_data_transfer/output_connector_data_transfer_linear_elasticity.tpp"
