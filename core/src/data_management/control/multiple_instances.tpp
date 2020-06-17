#include "data_management/control/multiple_instances.h"

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

template<typename FunctionSpaceType,typename BaseDataType>
MultipleInstances<FunctionSpaceType,BaseDataType>::
MultipleInstances(DihuContext context) : Data<FunctionSpaceType>(context)
{
}

template<typename FunctionSpaceType,typename BaseDataType>
MultipleInstances<FunctionSpaceType,BaseDataType>::
~MultipleInstances()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,typename BaseTimesteppingType>
void MultipleInstances<FunctionSpaceType,BaseTimesteppingType>::
setInstancesData(std::vector<BaseTimesteppingType> &instances)
{
  LOG(TRACE) << "setInstancesData for " << instances.size() << " instances";

  instancesData_.clear();
  instancesData_.reserve(instances.size());
  
  for (typename std::vector<BaseTimesteppingType>::iterator iter = instances.begin(); iter != instances.end(); iter++)
  {
    VLOG(1) << "MultipleInstancesData: push_back pointer to instance data";
    instancesData_.push_back(std::make_shared<BaseDataType>((*iter).data()));
  }

  // set functionSpace
  this->functionSpace_ = nullptr;
  if (!instances.empty())
  {
    VLOG(1) << " set Function space from instance data to MultipleInstances data";
    this->setFunctionSpace(instances[0].data().functionSpace());
  
    assert(this->functionSpace() != nullptr);
  }
}

template<typename FunctionSpaceType,typename BaseDataType>
void MultipleInstances<FunctionSpaceType,BaseDataType>::
createPetscObjects()
{
}

template<typename FunctionSpaceType,typename BaseDataType>
void MultipleInstances<FunctionSpaceType,BaseDataType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;
}

template<typename FunctionSpaceType,typename BaseDataType>
typename MultipleInstances<FunctionSpaceType,BaseDataType>::FieldVariablesForOutputWriter MultipleInstances<FunctionSpaceType,BaseDataType>::
getFieldVariablesForOutputWriter()
{
  std::vector<typename BaseDataType::FieldVariablesForOutputWriter> instancesFieldVariablesForOutputWriter;
  instancesFieldVariablesForOutputWriter.reserve(instancesData_.size());
  
  for (typename std::vector<std::shared_ptr<BaseDataType>>::const_iterator iter = instancesData_.begin(); iter != instancesData_.end(); iter++)
  {
    instancesFieldVariablesForOutputWriter.push_back((*iter)->getFieldVariablesForOutputWriter());
  }
  
  return std::make_tuple(instancesFieldVariablesForOutputWriter);
}


} // namespace
