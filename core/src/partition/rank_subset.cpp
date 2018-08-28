#include "partition/rank_subset.h"

#include <numeric>
#include <algorithm>

namespace Partition 
{
  
RankSubset::RankSubset()
{
  // create copy MPI_COMM_WORLD
  MPI_Comm_dup(MPI_COMM_WORLD, &mpiCommunicator_);
 
  // get number of ranks
  int nRanks;
  MPI_Comm_size(mpiCommunicator_, &nRanks);
  
  // create list of all ranks
  rankNo_.resize(nRanks);
  std::iota(rankNo_.begin(), rankNo_.end(), 0);
}
  
RankSubset::RankSubset(int singleRank)
{
  rankNo_.clear();
  rankNo_.push_back(singleRank);
  
  // get the own current MPI rank
  int currentRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);
  int color = MPI_UNDEFINED;
  
  // if currentRank is contained in rank subset
  if (singleRank == currentRank)
    color = 1;
  
  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_);
}

RankSubset::RankSubset(std::vector<int> &ranks) : rankNo_(ranks)
{
  // get the own current MPI rank
  int currentRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);
  int color = MPI_UNDEFINED;
  
  // if currentRank is contained in rank subset
  if (std::find(ranks.begin(),ranks.end(),currentRank) != ranks.end())
    color = 1;
  
  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_);
}

std::vector<int>::const_iterator RankSubset::begin()
{
  return rankNo_.cbegin();
}

std::vector<int>::const_iterator RankSubset::end()
{
  return rankNo_.cend();
}

element_no_t RankSubset::size() const
{
  return rankNo_.size();
}

MPI_Comm RankSubset::mpiCommunicator() const
{
  return mpiCommunicator_;
}
  
std::ostream &operator<<(std::ostream &stream, RankSubset rankSubset)
{
  if (rankSubset.size() == 0)
  {
    stream << "(empty rankSubset)";
  }
  else 
  {
    std::vector<int>::const_iterator iterRank = rankSubset.begin();
    stream << "(" << *iterRank;
    iterRank++;
    
    for (; iterRank != rankSubset.end(); iterRank++)
    {
      stream << ", " << *iterRank;
    }
    stream << ")";
  }
  return stream;
}

};
