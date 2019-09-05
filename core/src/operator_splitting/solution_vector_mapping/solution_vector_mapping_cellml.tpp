#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_cellml.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1>,
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>
>::transfer(const CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1> &transferableSolutionData1,
            const Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a>> fieldVariable1States        = transferableSolutionData1.stateVariable.values;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b>> fieldVariable1Intermediates = transferableSolutionData1.intermediateVariable.values;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>> fieldVariable2               = transferableSolutionData2.values;

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1States        = transferableSolutionData1.stateVariable.componentNo;
  int componentNo1Intermediates = transferableSolutionData1.intermediateVariable.componentNo;
  int componentNo2              = transferableSolutionData2.componentNo;

  double prefactor1 = transferableSolutionData1.stateVariable.scalingFactor;
  double prefactor2 = transferableSolutionData2.scalingFactor;

  bool useSlotIntermediates = (transferSlotName == "intermediates");

  // depending on the transferSlotName parameter, transfer states or intermediates
  if (useSlotIntermediates)
  {
    // componentNo1Intermediates is != -1, use the intermediates field variable
    VLOG(1) << "solution vector mapping (solution_vector_mapping_cellml.tpp), transfer intermediates from component "
      << componentNo1Intermediates << " (" << fieldVariable1Intermediates->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
      << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (not considered here)";

    assert(fieldVariable1Intermediates->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
    assert(fieldVariable1Intermediates->nDofsGlobal() == fieldVariable2->nDofsGlobal());

    // if representation of fieldVariable1Intermediates is invalid, this means that it has been extracted to another field variable
    if (fieldVariable2->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";
      // tranfer from finite elements back to cellml
      fieldVariable2->restoreExtractedComponent(fieldVariable1Intermediates->partitionedPetscVec());
    }
    else
    {
      // here we copy the given component of fieldVariable1Intermediates to the component of field variable2, because field variable 2 has > 1 components
      // if field fariable 2 has only 1 component, it could be extracted from fieldVariable1Intermediates without copy, this is done in the other specialization
      VLOG(1) << "SolutionVectorMapping VecCopy";
      PetscErrorCode ierr;
      ierr = VecCopy(fieldVariable1Intermediates->valuesGlobal(componentNo1Intermediates), fieldVariable2->valuesGlobal(componentNo2)); CHKERRV(ierr);
    }
  }
  else
  {
    // transfer the states field variable
    VLOG(1) << "solution vector mapping (solution_vector_mapping_cellml.tpp), transfer states from component "
      << componentNo1States << " (" << fieldVariable1States->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
      << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (not considered here)";

    assert(fieldVariable1States->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
    assert(fieldVariable1States->nDofsGlobal() == fieldVariable2->nDofsGlobal());

    // if representation of fieldVariable1States is invalid, this means that it has been extracted to another field variable
    if (fieldVariable2->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";
      // tranfer from finite elements back to cellml
      fieldVariable2->restoreExtractedComponent(fieldVariable1States->partitionedPetscVec());
    }
    else
    {
      // here we copy the given component of fieldVariable1States to the component of field variable2, because field variable 2 has > 1 components
      // if field fariable 2 has only 1 component, it could be extracted from fieldVariable1States without copy, this is done in the other specialization
      VLOG(1) << "SolutionVectorMapping VecCopy";
      PetscErrorCode ierr;
      ierr = VecCopy(fieldVariable1States->valuesGlobal(componentNo1States), fieldVariable2->valuesGlobal(componentNo2)); CHKERRV(ierr);
    }

    // scale result with prefactor1
    if (prefactor1 != 1.0)
    {
      PetscErrorCode ierr;
      ierr = VecScale(fieldVariable2->valuesGlobal(componentNo2), prefactor1); CHKERRV(ierr);
    }
  }
}

// reverse transfer
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void SolutionVectorMapping<
  Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>,
  CellMLOutputConnectorDataType<nComponents2a,nComponents2b,FunctionSpaceType2>
>::transfer(const Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
            const CellMLOutputConnectorDataType<nComponents2a,nComponents2b,FunctionSpaceType2> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>> fieldVariable1               = transferableSolutionData1.values;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a>> fieldVariable2States        = transferableSolutionData2.stateVariable.values;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b>> fieldVariable2Intermediates = transferableSolutionData2.intermediateVariable.values;

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1              = transferableSolutionData1.componentNo;
  int componentNo2States        = transferableSolutionData2.stateVariable.componentNo;
  int componentNo2Intermediates = transferableSolutionData2.intermediateVariable.componentNo;

  double prefactor1 = transferableSolutionData1.scalingFactor;
  double prefactor2 = transferableSolutionData2.stateVariable.scalingFactor;

  bool useSlotIntermediates = (transferSlotName == "intermediates");

  // depending on the transferSlotName parameter, transfer states or intermediates
  if (useSlotIntermediates)
  {
    // componentNo1 is != -1, use the intermediates field variable
    VLOG(1) << "solution vector mapping (solution_vector_mapping_cellml.tpp), transfer intermediates from component "
      << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
      << " to " << componentNo2Intermediates << " (" << fieldVariable2Intermediates->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (not considered here)";

    assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2Intermediates->nDofsLocalWithoutGhosts());
    assert(fieldVariable1->nDofsGlobal() == fieldVariable2Intermediates->nDofsGlobal());

    // if representation of fieldVariable1 is invalid, this means that it has been extracted to another field variable
    if (fieldVariable2Intermediates->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";
      // tranfer from finite elements back to cellml
      fieldVariable2Intermediates->restoreExtractedComponent(fieldVariable1->partitionedPetscVec());
    }
    else
    {
      // here we copy the given component of fieldVariable1 to the component of field variable2, because field variable 2 has > 1 components
      // if field fariable 2 has only 1 component, it could be extracted from fieldVariable1 without copy, this is done in the other specialization
      VLOG(1) << "SolutionVectorMapping VecCopy";
      PetscErrorCode ierr;
      ierr = VecCopy(fieldVariable1->valuesGlobal(componentNo1), fieldVariable2Intermediates->valuesGlobal(componentNo2Intermediates)); CHKERRV(ierr);
    }
  }
  else
  {
    // transfer the states field variable
    VLOG(1) << "solution vector mapping (solution_vector_mapping_cellml.tpp), transfer states from component "
      << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
      << " to " << componentNo2States << " (" << fieldVariable2States->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (not considered here)";

    assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2States->nDofsLocalWithoutGhosts());
    assert(fieldVariable1->nDofsGlobal() == fieldVariable2States->nDofsGlobal());

    // if representation of fieldVariable1 is invalid, this means that it has been extracted to another field variable
    if (fieldVariable2States->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";
      // transfer from finite elements back to cellml
      fieldVariable2States->restoreExtractedComponent(fieldVariable1->partitionedPetscVec());
    }
    else
    {
      // here we copy the given component of fieldVariable1 to the component of field variable2, because field variable 2 has > 1 components
      // if field fariable 2 has only 1 component, it could be extracted from fieldVariable1 without copy, this is done in the other specialization
      VLOG(1) << "SolutionVectorMapping VecCopy";
      PetscErrorCode ierr;
      ierr = VecCopy(fieldVariable1->valuesGlobal(componentNo1), fieldVariable2States->valuesGlobal(componentNo2States)); CHKERRV(ierr);
    }

    // scale result with prefactor1
    if (prefactor1 != 1.0)
    {
      PetscErrorCode ierr;
      ierr = VecScale(fieldVariable2States->valuesGlobal(componentNo2States), prefactor1); CHKERRV(ierr);
    }
  }
}

/** Transfer between two field variables, the first is vector-valued, use given component number, store in second, which is scalar
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2>
void SolutionVectorMapping<
  CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1>,
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,1>
>::transfer(const CellMLOutputConnectorDataType<nComponents1a,nComponents1b,FunctionSpaceType1> &transferableSolutionData1,
            const Data::ScaledFieldVariableComponent<FunctionSpaceType2,1> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a>> fieldVariable1States        = transferableSolutionData1.stateVariable.values;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b>> fieldVariable1Intermediates = transferableSolutionData1.intermediateVariable.values;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> fieldVariable2                          = transferableSolutionData2.values;

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1States        = transferableSolutionData1.stateVariable.componentNo;
  int componentNo1Intermediates = transferableSolutionData1.intermediateVariable.componentNo;
  int componentNo2 = 0;

  double prefactor1 = transferableSolutionData1.stateVariable.scalingFactor;
  double prefactor2 = transferableSolutionData2.scalingFactor;

  bool useSlotIntermediates = (transferSlotName == "intermediates");

  // depending on the transferSlotName parameter, transfer states or intermediates
  if (useSlotIntermediates)
  {
    // transfer the intermediates field variable
    VLOG(1) << "solution vector mapping, transfer from intermediates component "
      << componentNo1Intermediates << " (" << fieldVariable1Intermediates->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
      << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (2nd prefactor not considered here)";

    assert(fieldVariable1Intermediates->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
    assert(fieldVariable1Intermediates->nDofsGlobal() == fieldVariable2->nDofsGlobal());

    // if representation of fieldVariable1Intermediates is invalid, this means that it has been extracted to another field variable
    if (fieldVariable2->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";

      // transfer from finite elements back to cellml
      fieldVariable2->restoreExtractedComponent(fieldVariable1Intermediates->partitionedPetscVec());
    }
    else
    {
      VLOG(1) << "SolutionVectorMapping extractComponentShared";
      VLOG(2) << "original field variable: " << *fieldVariable1Intermediates;

      // fieldVariable2 has only 1 component
      // The following retrieves the raw memory pointer from the Petsc vector in fieldVariable1Intermediates and reuses it for fieldVariable2
      // that means that fieldVariable cannot be used anymore, only after restoreExtractedComponent was called on fieldVariable1Intermediates. This is done in the other solution_vector_mapping transfer call.
      fieldVariable1Intermediates->extractComponentShared(componentNo1Intermediates, fieldVariable2);


      VLOG(2) << "resulting field variable: " << *fieldVariable2;
    }
  }
  else
  {
    // transfer the states field variable
    VLOG(1) << "solution vector mapping, transfer from states component "
      << componentNo1States << " (" << fieldVariable1States->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
      << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (2nd prefactor not considered here)";

    assert(fieldVariable1States->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
    assert(fieldVariable1States->nDofsGlobal() == fieldVariable2->nDofsGlobal());

    // if representation of fieldVariable1States is invalid, this means that it has been extracted to another field variable
    if (fieldVariable2->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";

      // transfer from finite elements back to cellml
      fieldVariable2->restoreExtractedComponent(fieldVariable1States->partitionedPetscVec());
    }
    else
    {
      VLOG(1) << "SolutionVectorMapping extractComponentShared";
      VLOG(2) << "original field variable: " << *fieldVariable1States;

      // fieldVariable2 has only 1 component
      // The following retrieves the raw memory pointer from the Petsc vector in fieldVariable1States and reuses it for fieldVariable2
      // that means that fieldVariable cannot be used anymore, only after restoreExtractedComponent was called on fieldVariable1States. This is done in the other solution_vector_mapping transfer call.
      fieldVariable1States->extractComponentShared(componentNo1States, fieldVariable2);


      VLOG(2) << "resulting field variable: " << *fieldVariable2;
    }

    // scale result with prefactor1
    if (prefactor1 != 1.0)
    {
      PetscErrorCode ierr;
      ierr = VecScale(fieldVariable2->valuesGlobal(componentNo2), prefactor1); CHKERRV(ierr);
    }
  }
}

// reverse transfer not possible, because extracting can only be done to a 1-component field variable
