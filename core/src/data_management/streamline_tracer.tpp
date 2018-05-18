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

template<typename BasisOnMeshType,typename BaseDataType>
StreamlineTracer<BasisOnMeshType,BaseDataType>::
StreamlineTracer(DihuContext context) : Data<BasisOnMeshType>(context)
{
}

template<typename BasisOnMeshType,typename BaseDataType>
StreamlineTracer<BasisOnMeshType,BaseDataType>::
~StreamlineTracer()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
setBaseData(std::shared_ptr<BaseDataType> baseData)
{
  baseData_ = baseData;

  // set mesh
  this->setMesh(baseData_->mesh());
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
createPetscObjects()
{
  LOG(DEBUG)<<"StreamlineTracer<BasisOnMeshType,BaseDataType>::createPetscObjects()"<<std::endl;
  assert(this->mesh_);
  this->gradient_ = this->mesh_->template createFieldVariable<3>("gradient");
}

template<typename BasisOnMeshType,typename BaseDataType>
FieldVariable::FieldVariable<BasisOnMeshType,3> &StreamlineTracer<BasisOnMeshType,BaseDataType>::
gradient()
{
  return *this->gradient_;
}

template<typename BasisOnMeshType,typename BaseDataType>
dof_no_t StreamlineTracer<BasisOnMeshType,BaseDataType>::
nUnknowns()
{
  return this->mesh_->nNodes();
}

template<typename BasisOnMeshType,typename BaseDataType>
constexpr int StreamlineTracer<BasisOnMeshType,BaseDataType>::
getNDofsPerNode()
{
  return 1;
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4)<<"======================";

  int nEntries;
  VecGetSize(this->gradient_->values(), &nEntries);
  VLOG(4)<<"increment ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->gradient_->values());
  VLOG(4)<<"======================";
}

template<typename BasisOnMeshType,typename BaseDataType>
typename StreamlineTracer<BasisOnMeshType,BaseDataType>::OutputFieldVariables StreamlineTracer<BasisOnMeshType,BaseDataType>::
getOutputFieldVariables()
{
  return std::tuple_cat(baseData_->getOutputFieldVariables(),
                        std::tuple<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>>(gradient_));
}


} // namespace