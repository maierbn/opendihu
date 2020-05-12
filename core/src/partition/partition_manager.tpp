#include "partition/partition_manager.h"

#include <cstdlib>

#include "utility/mpi_utility.h"
#include "easylogging++.h"

namespace Partition
{

template<typename FunctionSpace>
std::shared_ptr<MeshPartition<FunctionSpace>> Manager::
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
  
  return std::make_shared<MeshPartition<FunctionSpace>>(nElementsGlobal, nNodesGlobal, nDofsGlobal, rankSubset);
}

// use nElementsLocal and nRanks, fill nElementsGlobal
template<typename FunctionSpace>
std::shared_ptr<MeshPartition<FunctionSpace>> Manager::
createPartitioningStructuredLocal(PythonConfig specificSettings, std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
                                  const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
                                  const std::array<int,FunctionSpace::dim()> nRanks, std::vector<int> rankNos)
{ 
  LOG(DEBUG) << "Partition::Manager::createPartitioningStructuredLocal from localSize " << nElementsLocal 
    << ", nRanks " << nRanks;
  
  const int D = FunctionSpace::dim();
  
  // the subset of ranks for the partition to be created
  std::shared_ptr<RankSubset> rankSubset;
  
  // if no nextRankSubset was specified, use all available ranks
  if (nextRankSubset_ == nullptr)
  {
    // create rank subset of all available MPI ranks
    rankSubset = std::make_shared<RankSubset>();
    LOG(DEBUG) << specificSettings.getStringPath() << ": create new rank subset " << *rankSubset;
  }
  else
  {
    // if nextRankSubset was specified, use it
    rankSubset = nextRankSubset_;
    LOG(DEBUG) << specificSettings.getStringPath() << ": use previously set rankSubset " << *rankSubset;
  }

  // if rank Nos to use were given, use them
  if (!rankNos.empty())
  {
    LOG(DEBUG) << specificSettings.getStringPath() << ": rankNos option was specified: " << rankNos;
    rankSubset = std::make_shared<RankSubset>(rankNos.begin(), rankNos.end(), rankSubset);
  }

  LOG(DEBUG) << specificSettings.getStringPath() << ": using rankSubset " << *rankSubset;

  if (!rankSubset->ownRankIsContained())
  {
    // return zero mesh partition
    nElementsGlobal.fill(0);
    std::array<global_no_t, D> beginGlobal({0});
    std::array<element_no_t,FunctionSpace::dim()> nElementsLocalZero({0});
    std::array<int,FunctionSpace::dim()> nRanksZero({0});
    return std::make_shared<MeshPartition<FunctionSpace>>(nElementsLocalZero, nElementsGlobal, beginGlobal, nRanksZero, rankSubset);
  }
  
  int rankNoSubsetCommunicator = rankSubset->ownRankNo();
  int nRanksSubsetCommunicator = rankSubset->size();

  int nRanksTotal = 1;
  for (int i = 0; i < D; i++)
  {
    nRanksTotal *= nRanks[i];
  }

  if (rankSubset->size() != nRanksTotal && !nextRankSubset_)
  {
    LOG(FATAL) << specificSettings.getStringPath() << ": mesh was created without `nextRankSubset_` being set by MultipleInstances. "
      << "Therefore, it is not known which ranks should be on this mesh. However, you could set \"rankNos\".";
  }
  
  if (nRanksSubsetCommunicator != nRanksTotal)
  {
    LOG(ERROR) << specificSettings.getStringPath() << ": Number of ranks (" << nRanksSubsetCommunicator << ") in rank subset does not match given nRanks in config " << nRanks << ", total " << nRanksTotal << ".";
  }
  
  std::array<int,3> rankGridCoordinate({0});  // the coordinate of the current rank in the nRanks[0] x nRanks[1] x nRanks[2] grid of ranks
  rankGridCoordinate[0] = rankNoSubsetCommunicator % nRanks[0];
  
  if (D >= 2) 
  {
    rankGridCoordinate[1] = int((rankNoSubsetCommunicator % (nRanks[0]*nRanks[1])) / nRanks[0]);
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
    VLOG(1) << "rankGridCoordinate: " << rankGridCoordinate << ", nRanks:" << nRanks << ", z color: " << oneDimensionCommunicatorColor;
  }
  
  // reduce the global sizes in the coordinate directions
  std::array<element_no_t,D> globalSizeMpi({0});   // note: because of MPI this cannot be of type global_no_t, but has to be the same as the send buffer
  std::array<element_no_t,D> ownGlobalSizeMpi = nElementsLocal;

  // add up number of local elements along x axis
  MPIUtility::handleReturnValue(MPI_Allreduce(&nElementsLocal[0], &ownGlobalSizeMpi[0], 1, MPIU_INT,
                                            MPI_SUM, oneDimensionCommunicator[0]));

  VLOG(1) << "reduce in x direction: " << nElementsLocal[0] << " -> " << ownGlobalSizeMpi[0];
  
  // add up number of local elements along y axis
  if (D >= 2)
  {
     MPIUtility::handleReturnValue(MPI_Allreduce(&nElementsLocal[1], &ownGlobalSizeMpi[1], 1, MPIU_INT,
                                               MPI_SUM, oneDimensionCommunicator[1]));

    VLOG(1) << "reduce in y direction: " << nElementsLocal[1] << " -> " << ownGlobalSizeMpi[1];
  }
  
  // add up number of local elements along z axis
  if (D >= 3)
  {
     MPIUtility::handleReturnValue(MPI_Allreduce(&nElementsLocal[2], &ownGlobalSizeMpi[2], 1, MPIU_INT,
                                               MPI_SUM, oneDimensionCommunicator[2]));

    VLOG(1) << "reduce in z direction: " << nElementsLocal[2] << " -> " << ownGlobalSizeMpi[2];
  }
  
  // now broadcast globalSizeMpi value to all ranks
  if (rankNoSubsetCommunicator == 0)
  {
    globalSizeMpi = ownGlobalSizeMpi;
  }

  // broadcast the global size from rank 0. This could be avoided, but it is done here to check if the value is the same on every rank.
  MPIUtility::handleReturnValue(MPI_Bcast(globalSizeMpi.data(), D, MPIU_INT, 0, rankSubset->mpiCommunicator()));

  bool globalSizeIsDifferent = ownGlobalSizeMpi[0] != globalSizeMpi[0];

  if (D >= 2)
    globalSizeIsDifferent = globalSizeIsDifferent || ownGlobalSizeMpi[1] != globalSizeMpi[1];

  if (D == 3)
    globalSizeIsDifferent = globalSizeIsDifferent || ownGlobalSizeMpi[2] != globalSizeMpi[2];

  if (globalSizeIsDifferent)
  {
    LOG(FATAL) << specificSettings.getStringPath() << ": The specified partitioning is invalid. "
      << "You set \"inputMeshIsGlobal\": False and specified "
      << "the local number of elements for each rank with \"nRanks\": " << nRanks << "\n"
      << "On rank " << rankNoSubsetCommunicator << " (rank grid coordinates " << rankGridCoordinate << ")"
      << ", the local number of elements is: " << nElementsLocal << ".\n"
      << "The determined global number of elements is: " << ownGlobalSizeMpi
      << ", however, rank 0 determines the global number of elements differently, as: " << globalSizeMpi << ".\n"
      << "Make sure that the \"nElements\" options on every rank define a regular grid in space.";
  }

  // compute beginGlobal values by prefix sum
  std::array<global_no_t, D> beginGlobal({0});
  std::array<global_no_t, D> nElementsLocalBigSize;
  for (int i = 0; i < D; i++)
  {
    nElementsLocalBigSize[i] = nElementsLocal[i];
  }

  MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocalBigSize[0], &beginGlobal[0], 1, MPI_LONG_LONG_INT, MPI_SUM, oneDimensionCommunicator[0]));
  
  if (D >= 2)
  {
    MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocalBigSize[1], &beginGlobal[1], 1, MPI_LONG_LONG_INT, MPI_SUM, oneDimensionCommunicator[1]));
  }
  
  if (D >= 3)
  {
    MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocalBigSize[2], &beginGlobal[2], 1, MPI_LONG_LONG_INT, MPI_SUM, oneDimensionCommunicator[2]));
  }
  
  for (int i = 0; i < D; i++)
  {
    nElementsGlobal[i] = globalSizeMpi[i];
    VLOG(1) << "set nElementsGlobal[" << i << "] = " << nElementsGlobal[i];
  }
  
  LOG(DEBUG) << "create new meshPartition, nElementsLocal: " << nElementsLocal << ", nElementsGlobal: " << nElementsGlobal
    << ", beginGlobal: " << beginGlobal << ", nRanks: " << nRanks << ", rankSubset : " << *rankSubset;

  // create a mesh partition with prescribed local partitions
  return std::make_shared<MeshPartition<FunctionSpace>>(nElementsLocal, nElementsGlobal, beginGlobal, nRanks, rankSubset);
}

// use nElementsGlobal, fill nElementsLocal and nRanks
template<typename FunctionSpace>
std::shared_ptr<MeshPartition<FunctionSpace>> Manager::
createPartitioningStructuredGlobal(PythonConfig specificSettings, const std::array<global_no_t,FunctionSpace::dim()> nElementsGlobal,
                                   std::array<element_no_t,FunctionSpace::dim()> &nElementsLocal,
                                   std::array<int,FunctionSpace::dim()> &nRanks, std::vector<int> rankNos)
{ 
  LOG(DEBUG) << "Partition::Manager::createPartitioningStructuredGlobal from nElementsGlobal " << nElementsGlobal;
  
  // the subset of ranks for the partition to be created
  std::shared_ptr<RankSubset> rankSubset;
  
  // if no nextRankSubset was specified, use all available ranks
  if (nextRankSubset_ == nullptr)
  {
    // create rank subset of all available MPI ranks
    rankSubset = std::make_shared<RankSubset>();
    LOG(DEBUG) << "why";
  }
  else 
  {
    // if nextRankSubset was specified, use it
    rankSubset = nextRankSubset_;
    LOG(DEBUG) << "this one better";
  }
  
  // if rank Nos to use were given, use them
  if (!rankNos.empty())
  {
    LOG(DEBUG) << "rankNos option was specified: " << rankNos;
    rankSubset = std::make_shared<RankSubset>(rankNos.begin(), rankNos.end(), rankSubset);
  }

  LOG(DEBUG) << "using rankSubset " << *rankSubset;
  
  // create meshPartition
  std::shared_ptr<MeshPartition<FunctionSpace>> meshPartition = std::make_shared<MeshPartition<FunctionSpace>>(nElementsGlobal, rankSubset);
  
  // set parameters localSize and nRanks
  for (int coordinateDirection = 0; coordinateDirection < FunctionSpace::dim(); coordinateDirection++)
  {
    nElementsLocal[coordinateDirection] = meshPartition->nElementsLocal(coordinateDirection);
    nRanks[coordinateDirection] = meshPartition->nRanks(coordinateDirection);
    VLOG(1) << "set nElementsLocal[" << coordinateDirection << "] = " << nElementsLocal[coordinateDirection];
    VLOG(1) << "set nRanks[" << coordinateDirection << "] = " << nRanks[coordinateDirection];
  }
  
  return meshPartition;
}

//! create new partitioning of a composite mesh, this emulates a normal mesh but the values are taken from the submeshes
template<typename BasisFunctionType, int D>
std::shared_ptr<MeshPartition<::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>>> Manager::
createPartitioningComposite(const std::vector<std::shared_ptr<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces, std::vector<int> rankNos)
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

  // if rank Nos to use were given, use them
  if (!rankNos.empty())
  {
    LOG(DEBUG) << "rankNos option was specified: " << rankNos;
    rankSubset = std::make_shared<RankSubset>(rankNos.begin(), rankNos.end(), rankSubset);
  }

  LOG(DEBUG) << "using rankSubset " << *rankSubset;

  typedef ::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType> FunctionSpace;

  // create meshPartition
  std::shared_ptr<MeshPartition<FunctionSpace>> meshPartition
    = std::make_shared<MeshPartition<FunctionSpace>>(subFunctionSpaces, rankSubset);

  return meshPartition;
}

}  // namespace
