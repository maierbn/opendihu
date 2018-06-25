#include "partition/rank_subset.h"

namespace Partition 
{
  
RankSubset::RankSubset(int singleRank)
{
  rankNo_.clear();
  rankNo_.push_back(singleRank);
  
  // get the own current MPI rank
  int currentRank;
  MPI_Comm_rank(MPI_COMM_WORLD, currentRank);
  int color = MPI_UNDEFINED;
  
  // if currentRank is contained in rank subset
  if (singleRank == currentRank)
    color = 1;
  
  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_);
}

RankSubset::RankSubset(std::vector<int> ranks) : rankNo_(ranks)
{
  // get the own current MPI rank
  int currentRank;
  MPI_Comm_rank(MPI_COMM_WORLD, currentRank);
  int color = MPI_UNDEFINED;
  
  // if currentRank is contained in rank subset
  if (ranks.find(currentRank) != ranks.end())
    color = 1;
  
  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_);
}

std::vector< int >::const_iterator RankSubset::begin()
{
  return rankNo_.cbegin();
}

std::vector< int >::const_iterator RankSubset::end()
{
  return rankNo_.cend();
}

element_no_t RankSubset::size()
{
  return rankNo_.size();
}

MPI_Comm RankSubset::mpiCommunicator()
{
  return mpiCommunicator_;
}
  
};