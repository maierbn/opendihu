#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "partition/01_mesh_partition.h"

namespace Partition
{

/** This manager object creates mesh partitions. They are needed for meshes/function spaces.
 *  The rank subset to use for the next partitioning to be generated can be specified by setRankSubsetForNextCreatedMesh.
  */
class Manager
{
public:
  
  //! constructor
  Manager(PyObject *settings);
 
  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh
  template<typename FunctionSpace>
  std::shared_ptr<MeshPartition<FunctionSpace>> createPartitioningUnstructured(global_no_t nElementsGlobal, global_no_t nNodesGlobal, global_no_t nDofsGlobal);

  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh, for a structured mesh, from global sizes
  //! use globalSize, fill localSize and nRanks
  template<typename FunctionSpace>
  std::shared_ptr<MeshPartition<FunctionSpace>> createPartitioningStructuredGlobal(const std::array<global_no_t,FunctionSpace::dim()> nElementsGlobal,
                                                                                 std::array<element_no_t,FunctionSpace::dim()> &nElementsLocal, 
                                                                                 std::array<int,FunctionSpace::dim()> &nRanks);

  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh, for a structured mesh, from local sizes
  //! use localSize and nRanks, fill globalSize
  //! @param nRanks The number of ranks in the coordinate directions.
  template<typename FunctionSpace>
  std::shared_ptr<MeshPartition<FunctionSpace>> createPartitioningStructuredLocal(std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
                                                                                const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
                                                                                const std::array<int,FunctionSpace::dim()> nRanks);

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