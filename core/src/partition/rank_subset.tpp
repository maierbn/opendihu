#include "partition/rank_subset.h"

#include <numeric>
#include <algorithm>

#include "utility/mpi_utility.h"
#include "utility/vector_operators.h"
#include "easylogging++.h"

namespace Partition 
{
  
template<typename Iter>
RankSubset::RankSubset(Iter ranksBegin, Iter ranksEnd) : ownRankNo_(-1)
{
  std::copy(ranksBegin, ranksEnd, std::inserter(rankNo_, rankNo_.begin()));

  // get the own rank in the MPI_WORLD communcator
  int ownRankWorldCommunicator;
  MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &ownRankWorldCommunicator), "MPI_Comm_rank");
  int color = MPI_UNDEFINED;
  
  // if ownRankWorldCommunicator is contained in rank subset
  if (std::find(ranksBegin,ranksEnd,ownRankWorldCommunicator) != ranksEnd)
    color = 1;
  
  VLOG(1) << "RankSubset constructor from rank list " << rankNo_ << ", ownRankWorldCommunicator=" << ownRankWorldCommunicator << ", color=" << color;

  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPIUtility::handleReturnValue(MPI_Comm_split(MPI_COMM_WORLD, color, 0, &mpiCommunicator_), "MPI_Comm_split");

  // update rankNo_ set
  if (color == 1)
  {
    int nRanksInCommunicator;
    MPIUtility::handleReturnValue(MPI_Comm_size(mpiCommunicator_, &nRanksInCommunicator));
    if (nRanksInCommunicator != rankNo_.size())
    {
      LOG(WARNING) << "Resizing nRanks from " << rankNo_.size() << " entr" << (rankNo_.size() == 1? "y" : "ies")
        << " to " << nRanksInCommunicator << " entr" << (nRanksInCommunicator == 1? "y" : "ies") << ", because the program is only run with so many ranks.";

      int i = 0;
      for (std::set<int>::iterator iter = rankNo_.begin(); iter != rankNo_.end(); iter++, i++)
      {
        if (i == nRanksInCommunicator)
        {
          // erase the rest of the entries
          rankNo_.erase(iter, rankNo_.end());
          break;
        }
      }
    }
  }

  // all ranks that are not part of the communcator will store "MPI_COMM_NULL" as mpiCommunicator_

  VLOG(1) << "RankSubset constructor for ranks " << rankNo_ << " done";
}

}  // namespace
