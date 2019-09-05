#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_fibers_emg.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

#include "mesh/mapping_between_meshes.h"
#include "mesh/mesh_manager.h"

template<typename BasisFunctionType, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents1a,nComponents1b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>,
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>
>::
transfer(const std::vector<std::vector<
           CellMLOutputConnectorDataType<
             nComponents1a,nComponents1b,
             FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
           >
         >> &transferableSolutionData1,
         Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2> transferableSolutionData2,
         const std::string transferSlotName)
{
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1a> FieldVariableType1a;
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1b> FieldVariableType1b;

  LOG(DEBUG) << "Solution vector mapping (solution_vector_mapping_fibers_emg.tpp)";

  // prepare the target mesh for the mapping, set all factors to zero
  typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2> FieldVariableType2;
  DihuContext::meshManager()->template prepareMapping<FieldVariableType2>(transferableSolutionData2.values);

  const int targetComponentNo = transferableSolutionData2.componentNo;

  bool useSlotIntermediates = (transferSlotName == "intermediates");

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData1.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData1[i].size(); j++)
    {
      std::shared_ptr<FieldVariableType1a> fieldVariable1States        = transferableSolutionData1[i][j].stateVariable.values;   // this is the transmembranePotential
      std::shared_ptr<FieldVariableType1b> fieldVariable1Intermediates = transferableSolutionData1[i][j].intermediateVariable.values;   // this is the stress

      const int sourceComponentNoStates        = transferableSolutionData1[i][j].stateVariable.componentNo;
      const int sourceComponentNoIntermediates = transferableSolutionData1[i][j].intermediateVariable.componentNo;

      if (useSlotIntermediates && sourceComponentNoIntermediates == -1)
      {
        LOG(FATAL) << "In transfer of OutputConnectorSlot, transferSlotName set to \"" << transferSlotName << "\" but sourceComponentNoIntermediates is -1.";
      }

      //LOG(DEBUG) << "transmembranePotential:" << *fieldVariable1;
      //LOG(FATAL) << "end";

      // map data between field variables

      // transfer states (Vm) or intermediates (gamma)
      if (useSlotIntermediates)
      {
        // transfer intermediates
        DihuContext::meshManager()->mapLowToHighDimension<FieldVariableType1b, FieldVariableType2>(
          fieldVariable1Intermediates, sourceComponentNoIntermediates, transferableSolutionData2.values, targetComponentNo);
      }
      else
      {
        // transfer state
        DihuContext::meshManager()->mapLowToHighDimension<FieldVariableType1a, FieldVariableType2>(
          fieldVariable1States, sourceComponentNoStates, transferableSolutionData2.values, targetComponentNo);
      }
    }
  }

  // finalize the mapping to the target mesh, compute final values by dividing by the factors
  DihuContext::meshManager()->template finalizeMapping<FieldVariableType2>(transferableSolutionData2.values, targetComponentNo);
}


template<typename FunctionSpaceType1, int nComponents1, typename BasisFunctionType, int nComponents2a, int nComponents2b>
void SolutionVectorMapping<
  Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>,  // <3D field variable>
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents2a,nComponents2b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>
>::
transfer(Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>,
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
