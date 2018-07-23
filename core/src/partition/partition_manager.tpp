#include "partition/partition_manager.h"

#include <cstdlib>

#include "utility/mpi_utility.h"
#include "easylogging++.h"

namespace Partition
{

template<typename BasisOnMesh>
std::shared_ptr<MeshPartition<BasisOnMesh>> Manager::
createPartitioning(global_no_t globalSize)
{ 
  // the subset of ranks for the partition to be created
  std::shared_ptr<RankSubset> rankSubset;
  
  // if no nextRankSubset was specified, use all available ranks
  if (nextRankSubset_ == nullptr)
  {
    // create rank subset of all available MPI ranks
    rankSubset = std::make_shared<RankSubset>();
  }
  else 
  {
    // if nextRankSubset was specified, use it
    rankSubset = nextRankSubset_;
  }
  
  return std::make_shared<MeshPartition<BasisOnMesh>>(globalSize, rankSubset);
}

template<int D, typename BasisOnMesh>
std::shared_ptr<MeshPartition<BasisOnMesh>> Manager::
createPartitioningStructured(std::array<global_no_t,D> globalSize)
{ 
  // the subset of ranks for the partition to be created
  std::shared_ptr<RankSubset> rankSubset;
  
  // if no nextRankSubset was specified, use all available ranks
  if (nextRankSubset_ == nullptr)
  {
    // create rank subset of all available MPI ranks
    rankSubset = std::make_shared<RankSubset>();
  }
  else 
  {
    // if nextRankSubset was specified, use it
    rankSubset = nextRankSubset_;
  }
  
  return std::make_shared<MeshPartition<BasisOnMesh>>(globalSize, rankSubset);
}

};    // namespace
