#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_finite_element_method_cellml.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

template<typename FunctionSpaceType1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,   // <fieldVariableType,componentNo>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>::transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
            const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double> &transferableSolutionData2)
{

  // call the transfer function of the parent class
  SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
  >::transfer(std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>(transferableSolutionData1, 0, 1.0),
              transferableSolutionData2);
  return;   // need to use the first implementation

  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> fieldVariable1 = &transferableSolutionData1;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>> fieldVariable2 = std::get<0>(transferableSolutionData2);

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  // tranfer from finite elements back to cellml
  fieldVariable1->restoreExtractedComponent(fieldVariable2->partitionedPetscVec());
}

template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2>
void SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>   // <fieldVariableType,componentNo>
>::transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double> &transferableSolutionData1,
            const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> &transferableSolutionData2)
{
  // call the transfer function of the parent class
  SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>, int, double>
  >::transfer(transferableSolutionData1,
              std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>, int, double>(transferableSolutionData2, 0, 1.0));
  return;

  // transfer from cellml to finite elements

  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>> fieldVariable1 = std::get<0>(transferableSolutionData1);
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> fieldVariable2 = transferableSolutionData2;

  // disable checking for nans and infs because it takes a lot of time
  //fieldVariable1->checkNansInfs();

  int componentNo1 = std::get<1>(transferableSolutionData1);
  double prefactor1 = std::get<2>(transferableSolutionData1);

  VLOG(1) << "solution vector mapping, transfer from " << fieldVariable1->name() << " component "
    << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
    << " to " << fieldVariable2->name();

  assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
  assert(fieldVariable1->nDofsGlobal() == fieldVariable2->nDofsGlobal());

  PetscErrorCode ierr;
  fieldVariable1->extractComponentShared(componentNo1, fieldVariable2);

  // scale result with prefactor1
  if (prefactor1 != 1.0)
  {
    ierr = VecScale(fieldVariable2->valuesGlobal(0), prefactor1); CHKERRV(ierr);
  }
}
