#include "output_connector_data_transfer/output_connector_data_transfer.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1>,   // <fieldVariableType,componentNo,prefactor>
  Data::ComponentOfFieldVariable<FunctionSpaceType2,nComponents2>
>::transfer(const Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
            Data::ComponentOfFieldVariable<FunctionSpaceType2,nComponents2> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1> FieldVariable1;
  typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2> FieldVariable2;

  std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1.values;
  std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2.values;

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1 = transferableSolutionData1.componentNo;
  int componentNo2 = transferableSolutionData2.componentNo;
  VLOG(1) << "solution vector mapping (output_connector_data_transfer.tpp), transfer from component "
    << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs)" <<
    << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs),";

  assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
  assert(fieldVariable1->nDofsGlobal() == fieldVariable2->nDofsGlobal());

  // perform the mapping
  DihuContext::meshManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
  DihuContext::meshManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, componentNo1, fieldVariable2, componentNo2);
  DihuContext::meshManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);

  /*
  // if representation of fieldVariable1 is invalid, this means that it has been extracted to another field variable
  if (fieldVariable2->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
  {
    VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";
    // tranfer from finite elements back to cellml
    fieldVariable2->restoreExtractedComponent(fieldVariable1->partitionedPetscVec());
  }
  else
  {
    // here we copy the given component of fieldVariable1 to the component of field variable2, because field variable 2 has > 1 components
    // if field fariable 2 has only 1 component, it could be extracted from fieldVariable1 without copy, this is done in the other specialization
    VLOG(1) << "SolutionVectorMapping VecCopy";
    PetscErrorCode ierr;
    ierr = VecCopy(fieldVariable1->valuesGlobal(componentNo1), fieldVariable2->valuesGlobal(componentNo2)); CHKERRV(ierr);
  }*/
}

/** Transfer between two field variables, the first is vector-valued, use given component number, store in second, which is scalar
 */

template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2>
void SolutionVectorMapping<
  Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1>,   // <fieldVariableType,componentNo,prefactor>
  Data::ComponentOfFieldVariable<FunctionSpaceType2,1>
>::transfer(const Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1> &transferableSolutionData1,
            Data::ComponentOfFieldVariable<FunctionSpaceType2,1> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1> FieldVariable1;
  typedef FieldVariable::FieldVariable<FunctionSpaceType2,1> FieldVariable2;

  std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1.values;
  std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2.values;

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1 = transferableSolutionData1.componentNo;
  int componentNo2 = 0;

  VLOG(1) << "solution vector mapping (output_connector_data_transfer.tpp), transfer from component "
    << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), " <<
    << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs),";

  assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
  assert(fieldVariable1->nDofsGlobal() == fieldVariable2->nDofsGlobal());

  // perform the mapping
  DihuContext::meshManager()->prepareMapping(fieldVariable1, fieldVariable2);
  DihuContext::meshManager()->map(fieldVariable1, componentNo1, fieldVariable2, componentNo2);
  DihuContext::meshManager()->finalizeMapping(fieldVariable1, fieldVariable2);
/*DihuContext::meshManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
  DihuContext::meshManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, componentNo1, fieldVariable2, componentNo2);
  DihuContext::meshManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);*/

  /*
  // if representation of fieldVariable1 is invalid, this means that it has been extracted to another field variable
  if (fieldVariable2->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
  {
    VLOG(1) << "SolutionVectorMapping restoreExtractedComponent";

    // transfer from finite elements back to cellml
    fieldVariable2->restoreExtractedComponent(fieldVariable1->partitionedPetscVec());
  }
  else
  {
    VLOG(1) << "SolutionVectorMapping extractComponentShared";
    VLOG(2) << "original field variable: " << *fieldVariable1;

    // fieldVariable2 has only 1 component
    // The following retrieves the raw memory pointer from the Petsc vector in fieldVariable1 and reuses it for fieldVariable2
    // that means that fieldVariable cannot be used anymore, only after restoreExtractedComponent was called on fieldVariable1. This is done in the other output_connector_data_transfer transfer call.
    fieldVariable1->extractComponentShared(componentNo1, fieldVariable2);

    VLOG(2) << "resulting field variable: " << *fieldVariable2;
  }*/
}
