#include "partition/partitioned_petsc_mat_base.h"


template<typename BasisOnMeshType>
PartitionedPetscMatBase<BasisOnMeshType>::
PartitionedPetscMatBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition, std::string name) :
  meshPartition_(meshPartition), name_(name)
{
}
  