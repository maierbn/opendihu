#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,   // <fieldVariableType,componentNo,prefactor>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>::transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double> &transferableSolutionData1,
            const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>> fieldVariable1 = std::get<0>(transferableSolutionData1);
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>> fieldVariable2 = std::get<0>(transferableSolutionData2);

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1 = std::get<1>(transferableSolutionData1);
  int componentNo2 = std::get<1>(transferableSolutionData2);

  double prefactor1 = std::get<2>(transferableSolutionData1);
  double prefactor2 = std::get<2>(transferableSolutionData2);

  VLOG(1) << "solution vector mapping (solution_vector_mapping.tpp), transfer from component "
    << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
    << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (not considered here)";

  assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
  assert(fieldVariable1->nDofsGlobal() == fieldVariable2->nDofsGlobal());

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
  }

  // scale result with prefactor1
  if (prefactor1 != 1.0)
  {
    PetscErrorCode ierr;
    ierr = VecScale(fieldVariable2->valuesGlobal(componentNo2), prefactor1); CHKERRV(ierr);
  }
}

/** Transfer between two field variables, the first is vector-valued, use given component number, store in second, which is scalar
 */

template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2>
void SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,   // <fieldVariableType,componentNo,prefactor>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>, int, double>
>::transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double> &transferableSolutionData1,
            const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>, int, double> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>> fieldVariable1 = std::get<0>(transferableSolutionData1);
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> fieldVariable2 = std::get<0>(transferableSolutionData2);

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1 = std::get<1>(transferableSolutionData1);
  int componentNo2 = 0;

  double prefactor1 = std::get<2>(transferableSolutionData1);
  double prefactor2 = std::get<2>(transferableSolutionData2);

  VLOG(1) << "solution vector mapping (solution_vector_mapping.tpp), transfer from component "
    << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
    << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (2nd prefactor not considered here)";

  assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
  assert(fieldVariable1->nDofsGlobal() == fieldVariable2->nDofsGlobal());

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
    // that means that fieldVariable cannot be used anymore, only after restoreExtractedComponent was called on fieldVariable1. This is done in the other solution_vector_mapping transfer call.
    fieldVariable1->extractComponentShared(componentNo1, fieldVariable2);


    VLOG(2) << "resulting field variable: " << *fieldVariable2;
  }

  // scale result with prefactor1
  if (prefactor1 != 1.0)
  {
    PetscErrorCode ierr;
    ierr = VecScale(fieldVariable2->valuesGlobal(componentNo2), prefactor1); CHKERRV(ierr);
  }
}
