#include "partition/partition_manager.h"

#include <cstdlib>

#include "utility/mpi_utility.h"
#include "partition/rank_subset.h"
#include "easylogging++.h"

namespace Partition
{

Manager::Manager(PyObject *specificSettings) : specificSettings_(specificSettings), nextRankSubset_(nullptr)
{
  MPIUtility::handleReturnValue (MPI_Comm_size(MPI_COMM_WORLD, &nRanksCommWorld_));
  MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNoCommWorld_));
}
  
int Manager::nRanksCommWorld()
{
  return nRanksCommWorld_;
}

int Manager::rankNoCommWorld()
{
  return rankNoCommWorld_;
}

void Manager::setRankSubsetForNextCreatedMesh(std::shared_ptr<RankSubset> nextRankSubset)
{
  nextRankSubset_ = nextRankSubset;
}

//! store the ranks which should be used for collective MPI operations
void Manager::setRankSubsetForCollectiveOperations(std::shared_ptr<RankSubset> rankSubset)
{
  rankSubsetForCollectiveOperations_ = rankSubset;
}

//! ranks which should be used for collective MPI operations
std::shared_ptr<RankSubset> Manager::rankSubsetForCollectiveOperations()
{
  if (rankSubsetForCollectiveOperations_ == nullptr)
  {
    // if rank subset was not set, constructs a rank subset with all ranks (MPICommWorld)
    rankSubsetForCollectiveOperations_ = std::make_shared<Partition::RankSubset>();
  }

  VLOG(1) << "get rankSubsetForCollectiveOperations " << *rankSubsetForCollectiveOperations_;

  return rankSubsetForCollectiveOperations_;
}

};    // namespace
