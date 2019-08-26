#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

/** Transfer between the output from multiple_instances and MultidomainSolver
 */
template<typename FieldVariableType1, typename FieldVariableType2>
class SolutionVectorMapping<
  std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>>,   // <fieldVariableType,componentNo,prefactor>
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>>  // <Petsc Vecs which are the sub-vectors of a nested vector, transmembranePotential>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>> &transferableSolutionData1,
                       std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>> transferableSolutionData2,
                       const std::string transferSlotName);
};

/** Transfer between the output from MultidomainSolver and MultipleInstances
 */
template<typename FieldVariableType1, typename FieldVariableType2>
class SolutionVectorMapping<
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>>,  // <Petsc Vecs which are the sub-vectors of a nested vector, transmembranePotential>
  std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>>   // <fieldVariableType,componentNo,prefactor>
>
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>> &transferableSolutionData2,
                       const std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>> &transferableSolutionData1,
                       const std::string transferSlotName);
};

#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_multidomain.tpp"
