#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_base.h"

template<typename FunctionSpaceType>
int PartitionedPetscVecBase<FunctionSpaceType>::vectorNo_ = 0;

template<typename FunctionSpaceType>
PartitionedPetscVecBase<FunctionSpaceType>::
PartitionedPetscVecBase(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, std::string name) :
  name_(name), meshPartition_(meshPartition), currentRepresentation_(Partition::values_representation_t::representationGlobal)
{
  // add a unique number to each vector in debug mode, to distinguish vectors with the same name
#ifndef NDEBUG
  std::stringstream nameStr;
  nameStr << name << "_" << vectorNo_;
  vectorNo_++;
  name_ = nameStr.str();
#endif
}

template<typename FunctionSpaceType>
std::vector<PetscInt> &PartitionedPetscVecBase<FunctionSpaceType>::
localDofNos()
{
  return this->meshPartition_->localDofNos();
}

//! get the meshPartition
template<typename FunctionSpaceType>
const std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> PartitionedPetscVecBase<FunctionSpaceType>::
meshPartition()
{
  return this->meshPartition_;
}

template<typename FunctionSpaceType>
std::string PartitionedPetscVecBase<FunctionSpaceType>::
name()
{
  return this->name_;
}

template<typename FunctionSpaceType>
Partition::values_representation_t PartitionedPetscVecBase<FunctionSpaceType>::
currentRepresentation()
{
  return currentRepresentation_;
}
