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

template<typename BasisOnMeshType,int nComponents>
TimeStepping<BasisOnMeshType,nComponents>::
TimeStepping(DihuContext context) : Data<BasisOnMeshType>(context)
{
}

template<typename BasisOnMeshType,int nComponents>
TimeStepping<BasisOnMeshType,nComponents>::
~TimeStepping()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
createPetscObjects()
{
  LOG(DEBUG)<<"TimeStepping<BasisOnMeshType,nComponents>::createPetscObjects("<<nComponents<<")"<<std::endl;
  assert(this->mesh_);
  this->solution_ = this->mesh_->template createFieldVariable<nComponents>("solution");
  this->increment_ = this->mesh_->template createFieldVariable<nComponents>("increment");
  //for the variant 1 of the implicit Euler
  //this->rhs_ = this->mesh_->template createFieldVariable<nComponents>("rhs");
}

template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeStepping<BasisOnMeshType,nComponents>::
solution()
{
  return *this->solution_;
}

template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeStepping<BasisOnMeshType,nComponents>::
increment()
{
  return *this->increment_;
}

/*
template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeStepping<BasisOnMeshType,nComponents>::
rhs()
{
  return *this->rhs_;
}
*/

template<typename BasisOnMeshType,int nComponents>
dof_no_t TimeStepping<BasisOnMeshType,nComponents>::
nUnknowns()
{
  return this->mesh_->nNodes() * nComponents;
}

template<typename BasisOnMeshType,int nComponents>
constexpr int TimeStepping<BasisOnMeshType,nComponents>::
getNDofsPerNode()
{
  return nComponents;
}

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4)<<"======================";

  int nEntries;
  VecGetSize(this->increment_->values(), &nEntries);
  VLOG(4)<<"increment ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->increment_->values());
  VLOG(4)<<"======================";

  VecGetSize(this->solution_->values(), &nEntries);
  VLOG(4)<<"solution ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->solution_->values());
  VLOG(4)<<"======================";
}

template<typename BasisOnMeshType,int nComponents>
typename TimeStepping<BasisOnMeshType,nComponents>::OutputFieldVariables TimeStepping<BasisOnMeshType,nComponents>::
getOutputFieldVariables()
{
  return OutputFieldVariables(
    std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->geometryField()),
    solution_
  );
}


} // namespace