#include "partition/partition_manager.h"

#include <cstdlib>

#include "utility/mpi_utility.h"
#include "easylogging++.h"

namespace Partition
{

Manager::Manager(PyObject *specificSettings) : specificSettings_(specificSettings), nextRankSubset_(nullptr)
{
  MPIUtility::handleReturnValue (MPI_Comm_size(MPI_COMM_WORLD, &nRanksCommWorld_));
  MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNoCommWorld_));
}
  
std::shared_ptr<Partition> Manager::createPartitioning(global_no_t globalSize)
{ 
  int nRanks;  // number of ranks in partition that will be created
  int rankNo;  // own rank no in partition that will be created 
  if (nextRankSubset_ == nullptr)
  {
    // no rank subset means use all ranks
    nRanks = nRanksCommWorld_;
    rankNo = rankNoCommWorld_;
  }
  if (nextRankSubset_ != nullptr)
  {
    nRanks = nextRankSubset_.size();
    MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNo));
  }
  
  global_no_t sizePerRank = globalSize / nRanks;
  global_no_t residual = globalSize - sizePerRank*nRanks;
  element_no_t localSize = element_no_t(sizePerRank);
  
  if (rankNo < residual)
    localSize++;
  
  global_no_t beginGlobal = rankNo * sizePerRank + std::min(residual, (global_no_t)rankNo);
  global_no_t endGlobal = beginGlobal + localSize;
  
  LOG(DEBUG) << "create partition global (" << beginGlobal << "," << endGlobal << ") with size local " << localSize << "/ global" << globalSize;
  
  return std::make_shared<Partition>(beginGlobal, endGlobal, localSize);
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

};    // namespace
