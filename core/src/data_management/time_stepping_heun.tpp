#include "data_management/time_stepping_heun.h"

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

template<typename BasisOnMeshType,int nComponents>
TimeSteppingHeun<BasisOnMeshType,nComponents>::
TimeSteppingHeun(DihuContext context) : TimeStepping<BasisOnMeshType,nComponents>::TimeStepping(context)
{
}

template<typename BasisOnMeshType,int nComponents>
TimeSteppingHeun<BasisOnMeshType,nComponents>::
~TimeSteppingHeun()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,int nComponents>
void TimeSteppingHeun<BasisOnMeshType,nComponents>::
createPetscObjects()
{
  TimeStepping<BasisOnMeshType,nComponents>::createPetscObjects();

  LOG(DEBUG) << "TimeSteppingHeun<BasisOnMeshType,nComponents>::createPetscObjects(" << nComponents << ")";
  
  this->intermediateIncrement_ = this->mesh_->template createFieldVariable<nComponents>("intermediateIncrement");
}

template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeSteppingHeun<BasisOnMeshType,nComponents>::
intermediateIncrement()
{
  return *this->intermediateIncrement_;
}

/*template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeSteppingHeun<BasisOnMeshType,nComponents>::
intermediateSolution()
{
  return *this->intermediateSolution_;
}*/


template<typename BasisOnMeshType,int nComponents>
void TimeSteppingHeun<BasisOnMeshType,nComponents>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->intermediateIncrement_;
  VLOG(4) << *this->increment_;
  VLOG(4) << *this->solution_;
  VLOG(4) << "======================";
}

} // namespace
