#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "partition/01_mesh_partition.h"

namespace Partition
{

class Manager
{
public:
  
  //! constructor
  Manager(PyObject *settings);
 
  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh
  template<typename BasisOnMesh>
  std::shared_ptr<MeshPartition<BasisOnMesh>> createPartitioning(global_no_t globalSize);

  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh, for a structured mesh, from global sizes
  //! use globalSize, fill localSize and nRanks
  template<typename BasisOnMesh>
  std::shared_ptr<MeshPartition<BasisOnMesh>> createPartitioningStructuredGlobal(const std::array<global_no_t,BasisOnMesh::dim()> globalSize, std::array<element_no_t,BasisOnMesh::dim()> &localSize, std::array<int,BasisOnMesh::dim()> &nRanks);

  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh, for a structured mesh, from local sizes
  //! use localSize and nRanks, fill globalSize
  //! @param nRanks The number of ranks in the coordinate directions.
  template<typename BasisOnMesh>
  std::shared_ptr<MeshPartition<BasisOnMesh>> createPartitioningStructuredLocal(std::array<global_no_t,BasisOnMesh::dim()> &globalSize, const std::array<element_no_t,BasisOnMesh::dim()> localSize, const std::array<int,BasisOnMesh::dim()> nRanks);

  //! store a rank subset that will be used for the next partitioning that will be created
  void setRankSubsetForNextCreatedMesh(std::shared_ptr<RankSubset> nextRankSubset);
  
  //! get the number of MPI ranks 
  int nRanksCommWorld();
  
  //! get the own MPI rank no
  int rankNoCommWorld();
  
private:
 
  PyObject *specificSettings_;  ///< the settings object for the partition manager
 
  int nRanksCommWorld_;  ///< number of MPI ranks, processes
  int rankNoCommWorld_;   ///< own process no.
  
  std::shared_ptr<RankSubset> nextRankSubset_;
};

};    // namespace

#include "partition_manager.tpp"