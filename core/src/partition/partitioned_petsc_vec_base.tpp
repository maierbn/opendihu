#include "partition/partitioned_petsc_vec_base.h"

template<typename BasisOnMeshType>
PartitionedPetscVecBase<BasisOnMeshType>::
PartitionedPetscVecBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition) :
  meshPartition_(meshPartition)
{
}

template<typename BasisOnMeshType>
std::vector<PetscInt> &PartitionedPetscVecBase<BasisOnMeshType>::
localNodeNos()
{
  return this->meshPartition_->localNodeNos();
}
