#include "partition/00_mesh_partition_base.h"

#include <algorithm>

namespace Partition
{
 
MeshPartitionBase::MeshPartitionBase(std::shared_ptr<RankSubset> rankSubset) :
  rankSubset_(rankSubset)
{
}

MeshPartitionBase::~MeshPartitionBase()
{
}

int MeshPartitionBase::nRanks() const
{
  return this->rankSubset_->size();
}

MPI_Comm MeshPartitionBase::mpiCommunicator() const
{
  return rankSubset_->mpiCommunicator();
}

}  // namespace
