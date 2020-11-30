#include "data_management/time_stepping/time_stepping_heun.h"

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

template<typename FunctionSpaceType,int nComponents>
TimeSteppingHeun<FunctionSpaceType,nComponents>::
TimeSteppingHeun(DihuContext context) : TimeStepping<FunctionSpaceType,nComponents>::TimeStepping(context)
{
  this->debuggingName_ = "Heun";
}

template<typename FunctionSpaceType,int nComponents>
TimeSteppingHeun<FunctionSpaceType,nComponents>::
~TimeSteppingHeun()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingHeun<FunctionSpaceType,nComponents>::
createPetscObjects()
{
  TimeStepping<FunctionSpaceType,nComponents>::createPetscObjects();

  LOG(DEBUG) << "TimeSteppingHeun<FunctionSpaceType,nComponents>::createPetscObjects(" << nComponents << ")";
  LOG(DEBUG) << "functionSpace has " << this->functionSpace_->nDofsLocalWithoutGhosts() << " local without ghost dofs";
  
  this->algebraicIncrement_ = this->functionSpace_->template createFieldVariable<nComponents>("algebraicIncrement");
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeSteppingHeun<FunctionSpaceType,nComponents>::
algebraicIncrement()
{
  return this->algebraicIncrement_;
}

/*template<typename FunctionSpaceType,int nComponents>
FieldVariable::FieldVariable<FunctionSpaceType,nComponents> &TimeSteppingHeun<FunctionSpaceType,nComponents>::
algebraicSolution()
{
  return *this->algebraicSolution_;
}*/


template<typename FunctionSpaceType,int nComponents>
void TimeSteppingHeun<FunctionSpaceType,nComponents>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->algebraicIncrement_;
  VLOG(4) << *this->increment_;
  VLOG(4) << *this->solution_;
  VLOG(4) << "======================";
}

} // namespace
