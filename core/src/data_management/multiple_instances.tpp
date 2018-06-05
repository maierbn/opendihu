#include "data_management/multiple_instances.h"

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
MultipleInstances<BasisOnMeshType,BaseDataType>::
MultipleInstances(DihuContext context) : Data<BasisOnMeshType>(context)
{
}

template<typename BasisOnMeshType,typename BaseDataType>
MultipleInstances<BasisOnMeshType,BaseDataType>::
~MultipleInstances()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename BaseTimesteppingType>
void MultipleInstances<BasisOnMeshType,BaseTimesteppingType>::
setInstancesData(std::vector<BaseTimesteppingType> &instances)
{
  instancesData_.clear();
  instancesData_.reserve(instances.size());
  
  for (typename std::vector<BaseTimesteppingType>::iterator iter = instances.begin(); iter != instances.end(); iter++)
  {
    instancesData_.push_back(std::make_shared<BaseDataType>((*iter).data()));
  }

  // set mesh
  this->setMesh(instances[0].data().mesh());
  
  assert(this->mesh() != nullptr);
}

template<typename BasisOnMeshType,typename BaseDataType>
void MultipleInstances<BasisOnMeshType,BaseDataType>::
createPetscObjects()
{
}

template<typename BasisOnMeshType,typename BaseDataType>
dof_no_t MultipleInstances<BasisOnMeshType,BaseDataType>::
nUnknowns()
{
  return this->mesh_->nNodes();
}

template<typename BasisOnMeshType,typename BaseDataType>
constexpr int MultipleInstances<BasisOnMeshType,BaseDataType>::
getNDofsPerNode()
{
  return 1;
}

template<typename BasisOnMeshType,typename BaseDataType>
void MultipleInstances<BasisOnMeshType,BaseDataType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;
}

template<typename BasisOnMeshType,typename BaseDataType>
typename MultipleInstances<BasisOnMeshType,BaseDataType>::OutputFieldVariables MultipleInstances<BasisOnMeshType,BaseDataType>::
getOutputFieldVariables()
{
  std::vector<typename BaseDataType::OutputFieldVariables> instancesOutputFieldVariables;
  instancesOutputFieldVariables.reserve(instancesData_.size());
  
  for (typename std::vector<std::shared_ptr<BaseDataType>>::const_iterator iter = instancesData_.begin(); iter != instancesData_.end(); iter++)
  {
    instancesOutputFieldVariables.push_back((*iter)->getOutputFieldVariables());
  }
  
  return std::make_tuple(instancesOutputFieldVariables);
}


} // namespace