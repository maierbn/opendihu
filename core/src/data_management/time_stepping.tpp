#include "data_management/time_stepping.h"

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

template<typename BasisOnMeshType>
TimeStepping<BasisOnMeshType>::
TimeStepping(const DihuContext &context) : Data<BasisOnMeshType>(context)
{
  this->disablePrinting_ = true;
  if (PythonUtility::containsKey(this->context_.getPythonConfig(), "disablePrinting"))
  {
    this->disablePrinting_ = PythonUtility::getOptionBool(this->context_.getPythonConfig(), "disablePrinting", false);
  }
}

template<typename BasisOnMeshType>
TimeStepping<BasisOnMeshType>::
~TimeStepping()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType>
void TimeStepping<BasisOnMeshType>::
createPetscObjects()
{
  LOG(DEBUG)<<"TimeStepping<BasisOnMeshType>::createPetscObjects("<<this->nComponentsPerNode_<<")"<<std::endl;
  this->solution_ = this->mesh_->createFieldVariable("solution", this->nComponentsPerNode_);
  this->increment_ = this->mesh_->createFieldVariable("increment", this->nComponentsPerNode_);
}

template<typename BasisOnMeshType>
FieldVariable::FieldVariable<BasisOnMeshType> &TimeStepping<BasisOnMeshType>::
solution()
{
  return *this->solution_;
}

template<typename BasisOnMeshType>
FieldVariable::FieldVariable<BasisOnMeshType> &TimeStepping<BasisOnMeshType>::
increment()
{
  return *this->increment_;
}

template<typename BasisOnMeshType>
void TimeStepping<BasisOnMeshType>::
print()
{
  if (this->disablePrinting_)
    return;
  
  VLOG(1)<<"======================";
  
  int nEntries;
  VecGetSize(this->increment_->values(), &nEntries);
  VLOG(1)<<"increment ("<<nEntries<<" entries):";
  VLOG(1)<<PetscUtility::getStringVector(this->increment_->values());
  VLOG(1)<<"======================";
  
  VecGetSize(this->solution_->values(), &nEntries);
  VLOG(1)<<"solution ("<<nEntries<<" entries):";
  VLOG(1)<<PetscUtility::getStringVector(this->solution_->values());
  VLOG(1)<<"======================";
}

template<typename BasisOnMeshType>
std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> TimeStepping<BasisOnMeshType>::
fieldVariables()
{
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> result;
  result.push_back(std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType>>(this->mesh_->geometryField()));
  result.push_back(solution_);
  result.push_back(increment_);
  
  return result;
}
  


} // namespace