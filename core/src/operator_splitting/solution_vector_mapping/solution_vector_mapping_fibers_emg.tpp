#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_fibers_emg.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

#include "mesh/mapping_between_meshes.h"
#include "mesh/mesh_manager.h"

template<typename BasisFunctionType, int nComponents1a, int nComponents1b, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents1a,nComponents1b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>,
  std::shared_ptr<FieldVariableType2>
>::
transfer(const std::vector<std::vector<
           CellMLOutputConnectorDataType<
             nComponents1a,nComponents1b,
             FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
           >
         >> &transferableSolutionData1,
         std::shared_ptr<FieldVariableType2> transferableSolutionData2,
         const std::string transferSlotName)
{
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1a> FieldVariableType1a;
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1b> FieldVariableType1b;

  LOG(DEBUG) << "Solution vector mapping (solution_vector_mapping_fibers_emg.tpp)";

  // prepare the target mesh for the mapping, set all factors to zero
  DihuContext::meshManager()->template prepareMapping<FieldVariableType2>(transferableSolutionData2);

  const int targetComponentNo = 0;

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData1.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData1[i].size(); j++)
    {
      std::shared_ptr<FieldVariableType1a> fieldVariable1States        = transferableSolutionData1[i][j].stateVariable.values;   // this is the transmembranePotential
      std::shared_ptr<FieldVariableType1b> fieldVariable1Intermediates = transferableSolutionData1[i][j].intermediateVariable.values;   // this is the stress

      const int sourceComponentNoStates        = transferableSolutionData1[i][j].stateVariable.componentNo;
      const int sourceComponentNoIntermediates = transferableSolutionData1[i][j].intermediateVariable.componentNo;

      //LOG(DEBUG) << "transmembranePotential:" << *fieldVariable1;
      //LOG(FATAL) << "end";

      // map data between field variables

      // depending on the value of sourceComponentNoIntermediates, transfer states (Vm) or intermediates (gamma)
      if (sourceComponentNoIntermediates == -1)
      {
        // transfer state
        DihuContext::meshManager()->map<FieldVariableType1a, FieldVariableType2>(fieldVariable1States, sourceComponentNoStates, transferableSolutionData2, targetComponentNo);
      }
      else
      {
        // transfer intermediates
        DihuContext::meshManager()->map<FieldVariableType1b, FieldVariableType2>(fieldVariable1Intermediates, sourceComponentNoIntermediates, transferableSolutionData2, targetComponentNo);
      }
    }
  }

  // finalize the mapping to the target mesh, compute final values by dividing by the factors
  DihuContext::meshManager()->template finalizeMapping<FieldVariableType2>(transferableSolutionData2, targetComponentNo);
}

template<typename BasisFunctionType, typename FieldVariableType1, int nComponents2a, int nComponents2b>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents2a,nComponents2b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>
>::
transfer(std::shared_ptr<FieldVariableType1> transferableSolutionData1,
         const std::vector<std::vector<
           CellMLOutputConnectorDataType<
             nComponents2a,nComponents2b,
             FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
           >
          >> &transferableSolutionData2,
         const std::string transferSlotName)
{
  // do nothing to map from 3D Vm field back to fibers
}
