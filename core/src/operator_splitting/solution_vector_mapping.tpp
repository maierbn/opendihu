#include "operator_splitting/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double>,   // <fieldVariableType,componentNo,prefactor>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>::transfer(std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>>, int, double> transferableSolutionData1,
            std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double> transferableSolutionData2)
{
  // rename input data
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1>> fieldVariable1 = std::get<0>(transferableSolutionData1);
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>> fieldVariable2 = std::get<0>(transferableSolutionData2);

  fieldVariable1->checkNansInfs();

  int componentNo1 = std::get<1>(transferableSolutionData1);
  int componentNo2 = std::get<1>(transferableSolutionData2);

  double prefactor1 = std::get<2>(transferableSolutionData1);
  double prefactor2 = std::get<2>(transferableSolutionData2);

  VLOG(1) << "solution vector mapping, transfer from component "
    << componentNo1 << " (" << fieldVariable1->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor1
    << " to " << componentNo2 << " (" << fieldVariable2->nDofsLocalWithoutGhosts() << " dofs), prefactor " << prefactor2 << " (not considered here)";

  assert(fieldVariable1->nDofsLocalWithoutGhosts() == fieldVariable2->nDofsLocalWithoutGhosts());
  assert(fieldVariable1->nDofsGlobal() == fieldVariable2->nDofsGlobal());

  PetscErrorCode ierr;
  ierr = VecCopy(fieldVariable1->valuesGlobal(componentNo1), fieldVariable2->valuesGlobal(componentNo2)); CHKERRV(ierr);

  // scale result with prefactor1
  if (prefactor1 != 1.0)
  {
    ierr = VecScale(fieldVariable2->valuesGlobal(componentNo2), prefactor1);
  }
}

template<typename FunctionSpaceType1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,   // <fieldVariableType,componentNo>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>::transfer(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> transferableSolutionData1,
            std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double> transferableSolutionData2)
{
  // call the transfer function of the parent class
  SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
  >::transfer(std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>(transferableSolutionData1, 0, 1.0),
              transferableSolutionData2);
}
