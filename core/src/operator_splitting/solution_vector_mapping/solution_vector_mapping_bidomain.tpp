#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

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
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType> FieldVariableType1;

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData1.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData1[i].size(); j++)
    {
       std::shared_ptr<FieldVariableType1> transmembranePotential = std::get<0>(transferableSolutionData1[i][j]);
       std::shared_ptr<MappingBetweenMeshes<FieldVariableType1::FunctionSpace, FieldVariableType2::FunctionSpace>> mappingBetweenMeshes
         = DihuContext::meshManager()->mappingBetweenMeshes(transmembranePotential->functionSpace()->meshName(), transferableSolutionData2->functionSpace()->meshName());
       mappingBetweenMeshes->map(transmembranePotential, 0, transferableSolutionData2, 0);
    }
  }
};

template<typename BasisFunctionType, typename FieldVariableType1>
void SolutionVectorMapping<
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1>>,
      int, double>
    >>,   // vector<vector<fieldVariableType,componentNo,prefactor>>
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>::transfer(std::shared_ptr<FieldVariableType2> transferableSolutionData1,
            std::vector<std::vector<
            std::tuple<
              std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, 1>>,
              int, double>
            >> &transferableSolutionData2)
{

};
