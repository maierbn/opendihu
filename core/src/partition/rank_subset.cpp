#include "partition/rank_subset.h"

#include <numeric>
#include <algorithm>

#include "utility/mpi_utility.h"
#include "utility/vector_operators.h"
#include "easylogging++.h"

namespace Partition 
{
  
RankSubset::RankSubset() : ownRankNo_(-1)
{
  VLOG(1) << "RankSubset empty constructor";

  // create copy MPI_COMM_WORLD
  MPIUtility::handleReturnValue(MPI_Comm_dup(MPI_COMM_WORLD, &mpiCommunicator_), "MPI_Comm_dup");
 
  // get number of ranks
  int nRanks;
  MPIUtility::handleReturnValue(MPI_Comm_size(mpiCommunicator_, &nRanks), "MPI_Comm_size");
  
  // create list of all ranks
  rankNo_.resize(nRanks);
  std::iota(rankNo_.begin(), rankNo_.end(), 0);
}
  
RankSubset::RankSubset(int singleRank) : ownRankNo_(-1)
{
  rankNo_.clear();
  rankNo_.push_back(singleRank);
  
  // get the own current MPI rank
  int currentRank;
  MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &currentRank), "MPI_Comm_rank");
  int color = MPI_UNDEFINED;
  
  // if currentRank is contained in rank subset
  if (singleRank == currentRank)
    color = 1;
  
  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPIUtility::handleReturnValue(MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_), "MPI_Comm_split");
}

RankSubset::RankSubset(std::vector<int> &ranks) : rankNo_(ranks), ownRankNo_(-1)
{
  // get the own current MPI rank
  int currentRank;
  MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &currentRank), "MPI_Comm_rank");
  int color = MPI_UNDEFINED;
  
  // if currentRank is contained in rank subset
  if (std::find(ranks.begin(),ranks.end(),currentRank) != ranks.end())
    color = 1;
  
  VLOG(1) << "RankSubset constructor from rank list " << rankNo_ << ", currentRank=" << currentRank << ", color=" << color;

  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPIUtility::handleReturnValue(MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_), "MPI_Comm_split");

  // update rankNo_ vector
  int nRanksInCommunicator;
  MPIUtility::handleReturnValue(MPI_Comm_size(mpiCommunicator_, &nRanksInCommunicator));
  if (nRanksInCommunicator != rankNo_.size())
  {
    LOG(WARNING) << "Resizing nRanks from " << rankNo_.size() << " entr" << (rankNo_.size() == 1? "y" : "ies")
      << " to " << nRanksInCommunicator << " entr" << (nRanksInCommunicator == 1? "y" : "ies") << ", because the program is only run with so many ranks.";
    rankNo_.resize(nRanksInCommunicator);
  }

  VLOG(1) << "RankSubset constructor done";
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

element_no_t RankSubset::ownRankNo()
{
  if (ownRankNo_ == -1)
  {
    // get the own rank id in this communicator
    MPIUtility::handleReturnValue(MPI_Comm_rank(mpiCommunicator_, &ownRankNo_), "MPI_Comm_rank");
  }
  return ownRankNo_;
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
