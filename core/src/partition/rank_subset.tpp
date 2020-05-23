#include "partition/rank_subset.h"

#include <numeric>
#include <algorithm>

#include "utility/mpi_utility.h"
#include "utility/vector_operators.h"
#include "easylogging++.h"

class DihuContext;   // forward declaration

namespace Partition 
{
  
template<typename Iter>
RankSubset::RankSubset(Iter ranksBegin, Iter ranksEnd, std::shared_ptr<RankSubset> parentRankSubset) : ownRankNo_(-1), nCommunicatorsSplit_(0)
{
  isWorldCommunicator_ = false;
  MPI_Comm parentCommunicator = MPI_COMM_WORLD;
  int ownRankParentCommunicator = 0;
  if (parentRankSubset)
  {
    LOG(DEBUG) << "create RankSubset with " << std::distance(ranksBegin,ranksEnd) << " ranks, use parent communicator";
    parentCommunicator = parentRankSubset->mpiCommunicator();
    ownRankParentCommunicator = parentRankSubset->ownRankNo();
    VLOG(1) << "determined own rank: " << ownRankParentCommunicator;
  }
  else
  {
    VLOG(1) << "create RankSubset with " << std::distance(ranksBegin,ranksEnd) << " ranks, from MPI_COMM_WORLD";
    // get the own rank in the parent communicator which is MPI_COMM_WORLD
    MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &ownRankParentCommunicator), "MPI_Comm_rank");
    VLOG(1) << "determined own rank: " << ownRankParentCommunicator;
  }

  // get the own rank in the communicator
  int color = MPI_UNDEFINED;
  
  std::copy(ranksBegin, ranksEnd, std::inserter(rankNo_, rankNo_.begin()));

  int orderKey = 0;
  // if ownRankParentCommunicator is contained in rank subset
  Iter pos = std::find(ranksBegin,ranksEnd,ownRankParentCommunicator);
  if (pos != ranksEnd)
  {
    color = 1;
    orderKey = std::distance(ranksBegin, pos);
  }

  VLOG(1) << "RankSubset constructor from rank list " << rankNo_ << ", ownRankParentCommunicator=" << ownRankParentCommunicator
    << ", color=" << color << ", orderKey: " << orderKey;

LOG(DEBUG) << "RankSubset constructor from rank list " << rankNo_ << ", ownRankParentCommunicator=" << ownRankParentCommunicator
    << ", color=" << color << ", orderKey: " << orderKey;

  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPIUtility::handleReturnValue(MPI_Comm_split(parentCommunicator, color, orderKey, &mpiCommunicator_), "MPI_Comm_split");

  // update rankNo_ set
  if (color == 1)
  {
    int nRanksInCommunicator;
    MPIUtility::handleReturnValue(MPI_Comm_size(mpiCommunicator_, &nRanksInCommunicator));
    LOG(DEBUG) << "nRanksInCommunicator" << nRanksInCommunicator;
    LOG(DEBUG) << "rankNo_.size()" << rankNo_.size();
    if (nRanksInCommunicator != rankNo_.size())
    {
      if (nRanksInCommunicator < rankNo_.size())
      {
        LOG(WARNING) << "Resizing nRanks from " << rankNo_.size() << " entr" << (rankNo_.size() == 1? "y" : "ies")
          << " to " << nRanksInCommunicator << " entr" << (nRanksInCommunicator == 1? "y" : "ies") << ", because the program is only run with so many ranks.";
      }
      else if (nRanksInCommunicator > rankNo_.size())
      {
        LOG(WARNING) << "Resizing nRanks from " << rankNo_.size() << " entr" << (rankNo_.size() == 1? "y" : "ies")
          << " to " << nRanksInCommunicator << " entr" << (nRanksInCommunicator == 1? "y" : "ies") << ", rankNo_: " << rankNo_ << ", parentRankSubset: ";

        if (parentRankSubset)
          LOG(INFO) << *parentRankSubset;
        else
          LOG(INFO) << "MPI_COMM_WORLD";

      }

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
#if 1
  // assign the name of the new communicator
  if (ownRankIsContained())
  {
    LOG(DEBUG) << "ownRankIsContained";

    // get name of old communicator
    std::vector<char> oldCommunicatorNameStr(MPI_MAX_OBJECT_NAME);
    int oldCommunicatorNameLength = 0;
    MPIUtility::handleReturnValue(MPI_Comm_get_name(parentCommunicator, oldCommunicatorNameStr.data(), &oldCommunicatorNameLength), "MPI_Comm_get_name");

    std::string oldCommunicatorName;
    if (oldCommunicatorNameLength > 0)
    {
      oldCommunicatorName = std::string(oldCommunicatorNameStr.begin(), oldCommunicatorNameStr.begin()+oldCommunicatorNameLength);
    }
    VLOG(1) << "oldCommunicatorName: " << oldCommunicatorName;

    // define new name
    std::stringstream newCommunicatorName;
    if (parentRankSubset)
    {
      newCommunicatorName << oldCommunicatorName << "_" << parentRankSubset->nCommunicatorsSplit();
    }
    else
    {
      newCommunicatorName << oldCommunicatorName << "_" << nWorldCommunicatorsSplit;
    }
    VLOG(1) << "newCommunicatorName: " << newCommunicatorName.str();

    // assign name
    communicatorName_ = newCommunicatorName.str();
    MPIUtility::handleReturnValue(MPI_Comm_set_name(mpiCommunicator_, communicatorName_.c_str()), "MPI_Comm_set_name");
  }
#endif
  if (parentRankSubset)
  {
    parentRankSubset->incrementNCommunicatorSplit();
  }
  else
  {
    nWorldCommunicatorsSplit++;
  }

  VLOG(1) << "RankSubset constructor for ranks " << rankNo_ << " done";
}

}  // namespace
