#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "cellml/00_cellml_adapter_base.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "output_connector_data_transfer/output_connection.h"

/** Transfer between the output from multiple_instances and MultidomainSolver
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FieldVariableType2>
class SolutionVectorMapping<
  std::vector<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>>,   // <fieldVariableType,componentNo,prefactor>
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>>  // <Petsc Vecs which are the sub-vectors of a nested vector, transmembranePotential>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::vector<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> &transferableSolutionData1,
                       std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>> transferableSolutionData2,
                       OutputConnection &outputConnection);
};

/** Transfer between the output from MultidomainSolver and MultipleInstances
 */
template<typename FieldVariableType1, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
class SolutionVectorMapping<
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>>,  // <Petsc Vecs which are the sub-vectors of a nested vector, transmembranePotential>
  std::vector<Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>>   // <fieldVariableType,componentNo,prefactor>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>> &transferableSolutionData1,
                       const std::vector<Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> &transferableSolutionData2,
                       OutputConnection &outputConnection);
};

#include "output_connector_data_transfer/output_connector_data_transfer_multidomain.tpp"
