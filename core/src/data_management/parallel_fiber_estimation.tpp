#include "data_management/parallel_fiber_estimation.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType>
ParallelFiberEstimation<FunctionSpaceType>::
ParallelFiberEstimation(DihuContext context) : Data<FunctionSpaceType>(context)
{
}

template<typename FunctionSpaceType>
ParallelFiberEstimation<FunctionSpaceType>::
~ParallelFiberEstimation()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType>
void ParallelFiberEstimation<FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "ParallelFiberEstimation<FunctionSpaceType>::createPetscObjects()" << std::endl;
  //assert(this->functionSpace_);
  
  // create field variables on local partition
  //this->gradient_ = this->functionSpace_->template createFieldVariable<3>("gradient");
}

template<typename FunctionSpaceType>
void ParallelFiberEstimation<FunctionSpaceType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  //VLOG(4) << *this->gradient_;
  VLOG(4) << "======================";
}

template<typename FunctionSpaceType>
typename ParallelFiberEstimation<FunctionSpaceType>::OutputFieldVariables ParallelFiberEstimation<FunctionSpaceType>::
getOutputFieldVariables()
{
  return std::tuple_cat(
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>>()
  );
}


} // namespace
