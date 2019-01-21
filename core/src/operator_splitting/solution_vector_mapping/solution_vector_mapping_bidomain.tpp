#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

#include "mesh/mapping_between_meshes.h"
#include "mesh/mesh_manager.h"

template<typename BasisFunctionType, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1>>,
      int, double>
    >>,   // vector<vector<fieldVariableType,componentNo,prefactor>>
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>:: transfer(const std::vector<std::vector<
             std::tuple<
               std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1>>,
               int, double>
             >> &transferableSolutionData1,
             std::shared_ptr<FieldVariableType2> transferableSolutionData2)
{
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1> FieldVariableType1;

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData1.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData1[i].size(); j++)
    {
      std::shared_ptr<FieldVariableType1> transmembranePotential = std::get<0>(transferableSolutionData1[i][j]);
      std::shared_ptr<Mesh::MappingBetweenMeshes<typename FieldVariableType1::FunctionSpace, typename FieldVariableType2::FunctionSpace>> mappingBetweenMeshes
        = std::static_pointer_cast<Mesh::MappingBetweenMeshes<typename FieldVariableType1::FunctionSpace, typename FieldVariableType2::FunctionSpace>>(
          DihuContext::meshManager()->mappingBetweenMeshes(transmembranePotential->functionSpace()->meshName(), transferableSolutionData2->functionSpace()->meshName())
        );
      if (!mappingBetweenMeshes)
      {

      }
      mappingBetweenMeshes->template map<1,FieldVariableType2::nComponents()>(*transmembranePotential, 0, *transferableSolutionData2, 0);
    }
  }
}

template<typename BasisFunctionType, typename FieldVariableType1>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1>>,
      int, double>
    >>   // vector<vector<fieldVariableType,componentNo,prefactor>>
>::transfer(std::shared_ptr<FieldVariableType1> transferableSolutionData1,
            std::vector<std::vector<
            std::tuple<
              std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1>>,
              int, double>
            >> transferableSolutionData2)
{
  // do nothing to map from 3D Vm field back to fibers
}
