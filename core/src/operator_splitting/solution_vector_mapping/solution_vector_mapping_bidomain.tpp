#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

#include "mesh/mapping_between_meshes/mapping_between_meshes.h"
#include "mesh/mesh_manager/mesh_manager.h"

template<typename BasisFunctionType, int nComponents1, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1>>,
      int, double>
    >>,   // vector<vector<fieldVariableType,componentNo,prefactor>>
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>:: transfer(const std::vector<std::vector<
             std::tuple<
               std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1>>,
               int, double>
             >> &transferableSolutionData1,
             std::shared_ptr<FieldVariableType2> transferableSolutionData2,
             const std::string transferSlotName)
{
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1> FieldVariableType1;

  VLOG(1) << " (solution_vector_mapping_bidomain.tpp)";

  // prepare the target mesh for the mapping, set all factors to zero
  DihuContext::meshManager()->template prepareMapping<FieldVariableType2>(transferableSolutionData2);

  const int targetComponentNo = 0;

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData1.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData1[i].size(); j++)
    {
      std::shared_ptr<FieldVariableType1> transmembranePotential = std::get<0>(transferableSolutionData1[i][j]);
      const int sourceComponentNo = std::get<1>(transferableSolutionData1[i][j]);

      //LOG(DEBUG) << "transmembranePotential:" << *transmembranePotential;
      //LOG(FATAL) << "end";

      // map data between field variables
      DihuContext::meshManager()->mapLowToHighDimension<FieldVariableType1, FieldVariableType2>(transmembranePotential, sourceComponentNo, transferableSolutionData2, targetComponentNo);
    }
  }

  // finalize the mapping to the target mesh, compute final values by dividing by the factors
  DihuContext::meshManager()->template finalizeMapping<FieldVariableType2>(transferableSolutionData2, targetComponentNo);
}

template<typename BasisFunctionType, typename FieldVariableType1, int nComponents2>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    std::tuple<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents2>>,
      int, double>
    >>   // vector<vector<fieldVariableType,componentNo,prefactor>>
>::transfer(std::shared_ptr<FieldVariableType1> transferableSolutionData1,
            std::vector<std::vector<
            std::tuple<
              std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents2>>,
              int, double>
            >> transferableSolutionData2,
            const std::string transferSlotName)
{
  // do nothing to map from 3D Vm field back to fibers
}
