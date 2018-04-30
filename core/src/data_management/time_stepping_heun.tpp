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
TimeSteppingHeun(DihuContext context) : TimeStepping(context)
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
  // initialize solution and increment variable by parental method
  TimeStepping<BasisOnMeshType,nComponents>::createPetscObjects();

  LOG(DEBUG)<<"TimeSteppingHeun<BasisOnMeshType,nComponents>::createPetscObjects("<<nComponents<<")"<<std::endl;
  this->intermediateIncrement_ = this->mesh_->template createFieldVariable<nComponents>("intermediateIncrement");
}

template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeSteppingHeun<BasisOnMeshType,nComponents>:: // hier unsicher ob "&TimeSteppingHeun" oder "&TimeStepping".
intermediateIncrement()
{
  return *this->intermediateIncrement_;
}


template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end. 
{
  if (!VLOG_IS_ON(4))
    return;
  
  VLOG(4)<<"======================";
  
  int nEntries;
  VecGetSize(this->intermediateIncrement_->values(), &nEntries);
  VLOG(4)<<"intermediateIncrement ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->intermediateIncrement_->values());
  VLOG(4)<<"======================";
  
  VecGetSize(this->increment_->values(), &nEntries);
  VLOG(4)<<"increment ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->increment_->values());
  VLOG(4)<<"======================";
  
  VecGetSize(this->solution_->values(), &nEntries);
  VLOG(4)<<"solution ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->solution_->values());
  VLOG(4)<<"======================";
}

} // namespace