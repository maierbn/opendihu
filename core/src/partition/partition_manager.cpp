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
