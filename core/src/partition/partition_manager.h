#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "partition/mesh_partition/01_mesh_partition.h"
#include "control/python_config/python_config.h"

// forward declaration
namespace FunctionSpace
{
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;
}

namespace Partition
{

/** This manager object creates mesh partitions. They are needed for meshes/function spaces.
 *  The rank subset to use for the next partitioning to be generated can be specified by setRankSubsetForNextCreatedPartitioning.
  */
class Manager
{
public:
  
  //! constructor
  Manager(PythonConfig settings);
 
  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedPartitioning
  template<typename FunctionSpace>
  std::shared_ptr<MeshPartition<FunctionSpace>> createPartitioningUnstructured(global_no_t nElementsGlobal, global_no_t nNodesGlobal, global_no_t nDofsGlobal);

  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedPartitioning, for a structured mesh, from global sizes
  //! use globalSize, fill localSize and nRanks
  template<typename FunctionSpace>
  std::shared_ptr<MeshPartition<FunctionSpace>> createPartitioningStructuredGlobal(const std::array<global_no_t,FunctionSpace::dim()> nElementsGlobal,
                                                                                 std::array<element_no_t,FunctionSpace::dim()> &nElementsLocal, 
                                                                                 std::array<int,FunctionSpace::dim()> &nRanks);

  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedPartitioning, for a structured mesh, from local sizes
  //! use localSize and nRanks, fill globalSize
  //! @param nRanks The number of ranks in the coordinate directions.
  template<typename FunctionSpace>
  std::shared_ptr<MeshPartition<FunctionSpace>> createPartitioningStructuredLocal(std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
                                                                                const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
                                                                                const std::array<int,FunctionSpace::dim()> nRanks);

  //! create new partitioning of a composite mesh, this emulates a normal mesh but the values are taken from the submeshes
  template<typename BasisFunctionType, int D>
  std::shared_ptr<MeshPartition<::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>>> createPartitioningComposite(
    const std::vector<std::shared_ptr<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces);

  //! store a rank subset that will be used for the next partitioning that will be created
  void setRankSubsetForNextCreatedPartitioning(std::shared_ptr<RankSubset> nextRankSubset);
  
  //! store the ranks which should be used for collective MPI operations
  void setRankSubsetForCollectiveOperations(std::shared_ptr<RankSubset> rankSubset);

  //! ranks which should be used for collective MPI operations
  std::shared_ptr<RankSubset> rankSubsetForCollectiveOperations();
  
private:
 
  PythonConfig specificSettings_;  ///< the settings object for the partition manager
  
  std::shared_ptr<RankSubset> nextRankSubset_;   ///< rank subset that will be used for the next partitioning that will be created
  std::shared_ptr<RankSubset> rankSubsetForCollectiveOperations_;    ///< the ranks which should be used for collective MPI operations
};

}  // namespace

#include "partition_manager.tpp"
