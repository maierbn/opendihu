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

void MeshPartitionBase::initializeLocalNodeNos(node_no_t nLocalDofsWithGhosts)
{
  // set the localDofNos_ vector to {0,1,2,...,nLocalDofsWithGhosts-1}
  localDofNos_.resize(nLocalDofsWithGhosts);
  std::iota(localDofNos_.begin(), localDofNos_.end(), 0);
}

int MeshPartitionBase::nRanks()
{
  return this->rankSubset_->size();
}

MPI_Comm MeshPartitionBase::mpiCommunicator()
{
  return rankSubset_->mpiCommunicator();
}

std::vector<PetscInt> &MeshPartitionBase::localDofNos()
{
  return localDofNos_;
}

}  // namespace
