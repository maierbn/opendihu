#include "output_connector_data_transfer/output_connector_data_transfer_linear_elasticity.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

#include "mesh/mapping_between_meshes/mapping_between_meshes.h"
#include "mesh/mesh_manager/mesh_manager.h"

template<typename BasisFunctionType, int nComponents1a, int nComponents1b, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents1a,nComponents1b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>,
  TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType2>
>::
transfer(const std::vector<std::vector<
           CellMLOutputConnectorDataType<
             nComponents1a,nComponents1b,
             FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
           >
         >> &transferableSolutionData1,
         TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType2> transferableSolutionData2,
         const std::string transferSlotName)
{
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1a> FieldVariableType1a;
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>, nComponents1b> FieldVariableType1b;

  LOG(DEBUG) << "Solution vector mapping (output_connector_data_transfer_linear_elasticity.tpp)";

  // prepare the target mesh for the mapping, set all factors to zero
  DihuContext::meshManager()->template prepareMapping<FieldVariableType2>(transferableSolutionData2.activation);

  const int targetComponentNo = 0; // 3D activation field variable has only 1 component

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
          fieldVariable1Intermediates, sourceComponentNoIntermediates, transferableSolutionData2.activation, targetComponentNo);
        LOG(DEBUG) << "transferred from         : " << *fieldVariable1Intermediates;
        LOG(DEBUG) << "transferred to activation: " << *transferableSolutionData2.activation;
      }
      else
      {
        // transfer state
        DihuContext::meshManager()->mapLowToHighDimension<FieldVariableType1a, FieldVariableType2>(
          fieldVariable1States, sourceComponentNoStates, transferableSolutionData2.activation, targetComponentNo);
      }
    }
  }

  // finalize the mapping to the target mesh, compute final values by dividing by the factors
  DihuContext::meshManager()->template finalizeMapping<FieldVariableType2>(transferableSolutionData2.activation, targetComponentNo);
}


template<typename FieldVariableType1, typename BasisFunctionType, int nComponents2a, int nComponents2b>
void SolutionVectorMapping<
  TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<
    CellMLOutputConnectorDataType<
      nComponents2a,nComponents2b,
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
    >
  >>
>::
transfer(TimeSteppingScheme::ElasticitySolverOutputConnectorDataType<FieldVariableType1> transferableSolutionData1,
         const std::vector<std::vector<
           CellMLOutputConnectorDataType<
             nComponents2a,nComponents2b,
             FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType>
           >
          >> &transferableSolutionData2,
         const std::string transferSlotName)
{
  // map geometry of 3D field back to all fiber geometries

  // define source and target field variable types
  typedef FieldVariable::FieldVariable<typename FieldVariableType1::FunctionSpace,3> SourceFieldVariable;

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunctionType> TargetFunctionSpace;
  typedef FieldVariable::FieldVariable<TargetFunctionSpace,3> TargetFieldVariable;

  // get source field variable, this is the same for all fibers
  std::shared_ptr<SourceFieldVariable> geometryFieldSource = std::make_shared<SourceFieldVariable>(transferableSolutionData1.activation->functionSpace()->geometryField());

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData2.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData2[i].size(); j++)
    {
      std::shared_ptr<TargetFieldVariable> geometryFieldTarget = std::make_shared<TargetFieldVariable>(transferableSolutionData2[i][j].stateVariable.values->functionSpace()->geometryField());

      // transfer values from 3D mesh to 1D mesh, using the inverse mapping that was already generated
      DihuContext::meshManager()->mapHighToLowDimension<SourceFieldVariable, TargetFieldVariable>(geometryFieldSource, geometryFieldTarget);
    }
  }
}
