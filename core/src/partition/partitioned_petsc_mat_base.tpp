#include "partition/partitioned_petsc_mat_base.h"


template<typename FunctionSpaceType>
PartitionedPetscMatBase<FunctionSpaceType>::
PartitionedPetscMatBase(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, std::string name) :
  meshPartition_(meshPartition), name_(name)
{
}
  