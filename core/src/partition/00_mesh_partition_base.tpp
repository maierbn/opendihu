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

void MeshPartitionBase::initializeLocalDofs()
{
  // set the localDofs_ vector to {0,1,2,...,localSize()-1}
  localDofs_.resize(localSize());
  std::iota(localDofs_.begin(), localDofs_.end(), 0);
}

int MeshPartitionBase::nRanks()
{
  return this->rankSubset_->size();
}

MPI_Comm MeshPartitionBase::mpiCommunicator()
{
  return rankSubset_->mpiCommunicator();
}

std::vector<PetscInt> &MeshPartitionBase::localDofs()
{
  return localDofs_;
}

}  // namespace