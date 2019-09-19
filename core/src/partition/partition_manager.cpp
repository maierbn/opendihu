#include "partition/partition_manager.h"

#include <cstdlib>

#include "utility/mpi_utility.h"
#include "partition/rank_subset.h"
#include "easylogging++.h"

namespace Partition
{

Manager::Manager(PythonConfig specificSettings) : specificSettings_(specificSettings), nextRankSubset_(nullptr)
{
}

void Manager::setRankSubsetForNextCreatedPartitioning(std::shared_ptr<RankSubset> nextRankSubset)
{
  nextRankSubset_ = nextRankSubset;
}

//! store the ranks which should be used for collective MPI operations
void Manager::setRankSubsetForCollectiveOperations(std::shared_ptr<RankSubset> rankSubset)
{
  // if the rank subset for collective operations was not yet set, set it now
  if (!rankSubsetForCollectiveOperations_)
  {
    rankSubsetForCollectiveOperations_ = rankSubset;
  }
  else if (rankSubsetForCollectiveOperations_->size() < rankSubset->size())
  {
    // if the rank subset for collective operations was set earlier, only overwrite if with a rank subset that includes more ranks
    // this happens when more MultipleInstances are nested
    rankSubsetForCollectiveOperations_ = rankSubset;
  }
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

}  // namespace
