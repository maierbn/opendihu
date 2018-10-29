#include "data_management/solution_vector_mapping.h"

#include <vector>
#include "easylogging++.h"

template<typename FieldVariable1, typename FieldVariable2>
void SolutionVectorMapping::
transfer(FieldVariable1 &solution1, std::shared_ptr<SolutionVectorMapping> solutionVectorMapping2, FieldVariable2 &solution2)
{
  VLOG(1) << "solution vector mapping, transfer from component "
    << outputComponentNo_ << " (" << solution1.nDofsLocalWithoutGhosts() << " dofs)"
    << " to " << solutionVectorMapping2->outputComponentNo_ << " (" << solution2.nDofsLocalWithoutGhosts() << " dofs)";

  assert(solution1.nDofsLocalWithoutGhosts() == solution2.nDofsLocalWithoutGhosts());
  assert(solution1.nDofsGlobal() == solution2.nDofsGlobal());

  PetscErrorCode ierr;
  ierr = VecCopy(solution1.valuesGlobal(outputComponentNo_), solution2.valuesGlobal(solutionVectorMapping2->outputComponentNo_)); CHKERRV(ierr);

  // scale data with scalingFactor
  if (scalingFactor_ != 1.0)
  {
    VecScale(solution2.valuesGlobal(solutionVectorMapping2->outputComponentNo_), scalingFactor_);
  }

  solution2.checkNansInfs();
}
