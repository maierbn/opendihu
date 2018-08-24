#include "partition/partitioned_petsc_vec_base.h"

template<typename BasisOnMeshType>
PartitionedPetscVecBase<BasisOnMeshType>::
PartitionedPetscVecBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition, std::string name) :
  name_(name), meshPartition_(meshPartition)
{
}

template<typename BasisOnMeshType>
std::vector<PetscInt> &PartitionedPetscVecBase<BasisOnMeshType>::
localDofNos()
{
  return this->meshPartition_->localDofNos();
}

//! get the meshPartition
template<typename BasisOnMeshType>
const std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> PartitionedPetscVecBase<BasisOnMeshType>::
meshPartition()
{
  return this->meshPartition_;
}

template<typename BasisOnMeshType>
std::string PartitionedPetscVecBase<BasisOnMeshType>::
name()
{
  return this->name_;
}

