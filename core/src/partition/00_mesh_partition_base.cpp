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

void MeshPartitionBase::initializeLocalNodeNos(node_no_t nLocalNodes)
{
  // set the localNodeNos_ vector to {0,1,2,...,nLocalNodes-1}
  localNodeNos_.resize(localSize());
  std::iota(localNodeNos_.begin(), localNodeNos_.end(), 0);
}

int MeshPartitionBase::nRanks()
{
  return this->rankSubset_->size();
}

MPI_Comm MeshPartitionBase::mpiCommunicator()
{
  return rankSubset_->mpiCommunicator();
}

std::vector<PetscInt> &MeshPartitionBase::localNodeNos()
{
  return localNodeNos_;
}

}  // namespace
