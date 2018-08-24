#include "partition/partition_manager.h"

#include <cstdlib>

#include "utility/mpi_utility.h"
#include "easylogging++.h"

namespace Partition
{

template<typename BasisOnMesh>
std::shared_ptr<MeshPartition<BasisOnMesh>> Manager::
createPartitioningUnstructured(global_no_t nElementsGlobal, global_no_t nNodesGlobal, global_no_t nDofsGlobal)
{ 
  LOG(DEBUG) << "Partition::Manager::createPartitioningUnstructured, nElementsGlobal: " 
    << nElementsGlobal << ", nNodesGlobal: " << nNodesGlobal << ", nDofsGlobal: " << nDofsGlobal;
  
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
  
  LOG(DEBUG) << "using rankSubset " << *rankSubset;
  
  return std::make_shared<MeshPartition<BasisOnMesh>>(nElementsGlobal, nNodesGlobal, nDofsGlobal, rankSubset);
}

// use nElementsLocal and nRanks, fill nElementsGlobal
template<typename BasisOnMesh>
std::shared_ptr<MeshPartition<BasisOnMesh>> Manager::
createPartitioningStructuredLocal(std::array<global_no_t,BasisOnMesh::dim()> &nElementsGlobal,
                                  const std::array<element_no_t,BasisOnMesh::dim()> nElementsLocal,
                                  const std::array<int,BasisOnMesh::dim()> nRanks)
{ 
  LOG(DEBUG) << "Partition::Manager::createPartitioningStructuredLocal from localSize " << nElementsLocal 
    << ", nRanks " << nRanks;
  
  const int D = BasisOnMesh::dim();
  
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
  
  int rankNoSubsetCommunicator;
  int nRanksSubsetCommunicator;
  MPIUtility::handleReturnValue(MPI_Comm_rank(rankSubset->mpiCommunicator(), &rankNoSubsetCommunicator));
  MPIUtility::handleReturnValue(MPI_Comm_size(rankSubset->mpiCommunicator(), &nRanksSubsetCommunicator));
  int nRanksTotal = 0;
  for (int i = 0; i < D; i++)
  {
    nRanksTotal += nRanks[i];
  }
  
  if (nRanksSubsetCommunicator != nRanksTotal)
  {
    LOG(ERROR) << "Number of ranks (" << nRanksSubsetCommunicator << ") in rank subset does not match given nRanks in config " << nRanks << " = " << nRanksTotal << ".";
  }
  
  std::array<int,3> rankGridCoordinate({0});  // the coordinate of the current rank in the nRanks[0] x nRanks[1] x nRanks[2] grid of ranks
  rankGridCoordinate[0] = rankNoSubsetCommunicator % nRanks[0];
  
  if (D >= 2) 
  {
    rankGridCoordinate[1] = int(rankNoSubsetCommunicator / nRanks[0]);
  }
  if (D >= 3)
  {
    rankGridCoordinate[2] = int(rankNoSubsetCommunicator / (nRanks[0]*nRanks[1]));
  }
  
  // expand nRanks to 3 entries where not valid entries are set to 1
  std::array<int,3> nRanks3({1});
  for (int i = 0; i < D; i++)
  {
    nRanks3[i] = nRanks[i];
  }
  
  // create communicators that contain only one row of elements
  std::array<MPI_Comm,D> oneDimensionCommunicator;  // 3 communicators (for D==3) that contain one row of elements for x,y and z coordinate directions
  
  // communicator of row in x-direction
  int oneDimensionCommunicatorColor = rankGridCoordinate[2]*nRanks[1] + rankGridCoordinate[1];
    
  // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
  MPIUtility::handleReturnValue(MPI_Comm_split(rankSubset->mpiCommunicator(), oneDimensionCommunicatorColor, rankNoSubsetCommunicator,
                 &oneDimensionCommunicator[0]));
  
  if (D >= 2)
  {
    // communicator of row in y-direction
    oneDimensionCommunicatorColor = rankGridCoordinate[2]*nRanks[0] + rankGridCoordinate[0];
   
    // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
    MPIUtility::handleReturnValue(MPI_Comm_split(rankSubset->mpiCommunicator(), oneDimensionCommunicatorColor, rankNoSubsetCommunicator,
                   &oneDimensionCommunicator[1]));
  }
  
  if (D >= 3)
  {
    // communicator of row in z-direction
    oneDimensionCommunicatorColor = rankGridCoordinate[1]*nRanks[0] + rankGridCoordinate[0];
   
    // create new communicator which contains all ranks that have the same value of color (and not MPI_UNDEFINED)
    MPIUtility::handleReturnValue(MPI_Comm_split(rankSubset->mpiCommunicator(), oneDimensionCommunicatorColor, rankNoSubsetCommunicator,
                   &oneDimensionCommunicator[2]));
  }
  
  // reduce the global sizes in the coordinate directions
  std::array<element_no_t,D> globalSizeMpi;   // note: because of MPI this cannot be of type global_no_t, but has to be the same as the send buffer
  if (rankGridCoordinate[1] == 0 && rankGridCoordinate[2] == 0)
  {
    MPIUtility::handleReturnValue(MPI_Reduce(&nElementsLocal[0], &globalSizeMpi[0], 1, MPI_INT, 
                                              MPI_SUM, 0, oneDimensionCommunicator[0]));
  }
  
  if (D >= 2 && rankGridCoordinate[0] == 0 && rankGridCoordinate[2] == 0)
  {
     MPIUtility::handleReturnValue(MPI_Reduce(&nElementsLocal[1], &globalSizeMpi[1], 1, MPI_INT, 
                                               MPI_SUM, 0, oneDimensionCommunicator[1]));
  }
  
  if (D >= 3 && rankGridCoordinate[0] == 0 && rankGridCoordinate[1] == 0)
  {
     MPIUtility::handleReturnValue(MPI_Reduce(&nElementsLocal[2], &globalSizeMpi[2], 1, MPI_INT, 
                                               MPI_SUM, 0, oneDimensionCommunicator[2]));
  }
  
  // now broadcast globalSizeMpi value to all ranks
  MPIUtility::handleReturnValue(MPI_Bcast(globalSizeMpi.data(), D, MPI_INT, 0, rankSubset->mpiCommunicator()));
  
  // compute beginGlobal values by prefix sum
  std::array<PetscInt, D> beginGlobal({0});
  MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocal[0], &beginGlobal[0], 1, MPI_INT, MPI_SUM, oneDimensionCommunicator[0]));
  
  if (D >= 2)
  {
    MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocal[1], &beginGlobal[1], 1, MPI_INT, MPI_SUM, oneDimensionCommunicator[1]));
  }
  
  if (D >= 3)
  {
    MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocal[2], &beginGlobal[2], 1, MPI_INT, MPI_SUM, oneDimensionCommunicator[2]));
  }
  
  for (int i = 0; i < D; i++)
  {
    nElementsGlobal[i] = globalSizeMpi[i];
    VLOG(1) << "set nElementsGlobal[" << i << "] = " << nElementsGlobal[i];
  }
  
  // create a mesh partition with prescribed local partitions
  return std::make_shared<MeshPartition<BasisOnMesh>>(nElementsLocal, nElementsGlobal, beginGlobal, nRanks, rankSubset);
}

// use nElementsGlobal, fill nElementsLocal and nRanks
template<typename BasisOnMesh>
std::shared_ptr<MeshPartition<BasisOnMesh>> Manager::
createPartitioningStructuredGlobal(const std::array<global_no_t,BasisOnMesh::dim()> nElementsGlobal, 
                                   std::array<element_no_t,BasisOnMesh::dim()> &nElementsLocal,
                                   std::array<int,BasisOnMesh::dim()> &nRanks)
{ 
  LOG(DEBUG) << "Partition::Manager::createPartitioningStructuredGlobal from nElementsGlobal " << nElementsGlobal;
  
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
  
  LOG(DEBUG) << "using rankSubset " << *rankSubset;
  
  // create meshPartition
  std::shared_ptr<MeshPartition<BasisOnMesh>> meshPartition = std::make_shared<MeshPartition<BasisOnMesh>>(nElementsGlobal, rankSubset);
  
  // set parameters localSize and nRanks
  for (int coordinateDirection = 0; coordinateDirection < BasisOnMesh::dim(); coordinateDirection++)
  {
    nElementsLocal[coordinateDirection] = meshPartition->nElementsLocal(coordinateDirection);
    nRanks[coordinateDirection] = meshPartition->nRanks(coordinateDirection);
    VLOG(1) << "set nElementsLocal[" << coordinateDirection << "] = " << nElementsLocal[coordinateDirection];
    VLOG(1) << "set nRanks[" << coordinateDirection << "] = " << nRanks[coordinateDirection];
  }
  
  return meshPartition;
}

};    // namespace
