#pragma once

#include <memory>
#include <petscdmda.h>

#include "partition/00_mesh_partition_base.h"
#include "control/types.h"
#include "partition/rank_subset.h"


namespace Partition
{
 
/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: ghost elements
 */
template<typename BasisOnMeshType, typename = typename BasisOnMeshType::Mesh>
class MeshPartition
{
};

/** partial specialization for structured meshes */
template<int D, typename MeshType, typename BasisFunctionType>
class MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>> :
  public MeshPartitionBase
{
public:
 
  //! constructor, determine the decomposition by PETSc
  MeshPartition(std::array<node_no_t,D> globalSize, std::shared_ptr<RankSubset> rankSubset);
 
  //! constructor from prescribed partition
  MeshPartition(std::array<node_no_t,D> localSize, std::array<node_no_t,D> globalSize, 
                std::array<int,D> nRanks, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of ranks in a coordinate direction
  int nRanks(int coordinateDirection);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t localSize();
  
  //! number of entries in the given coordinate direction in the current partition
  element_no_t localSize(int coordinateDirection);

  //! first local number in current partition
  global_no_t beginGlobal(int coordinateDirection);
  
  //! one after last number in current partition
  global_no_t endGlobal(int coordinateDirection);
  
  //! number of nodes in total
  global_no_t globalSize();
  
  //! number of nodes in total, in the given coordinate direction 
  global_no_t globalSize(int coordinateDirection);
  
  //! get if the local partition touches the right/top/back border
  bool localPartitionIsAtBorder(coordinateDirection);
  
  //! get a vector with the local sizes on every rank
  std::vector<element_no_t> &localSizesOnRanks(int coordinateDirection);
  
  //! get an AO object TODO: is this needed anywhere? If not, remove
  AO &applicationOrdering();
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! from a vector of global numbers remove all that are non-local
  template <typename T>
  void extractLocalNumbers(std::vector<T> &vector);
  
protected:
 
  DM dm_;    ///< PETSc DMDA object (data management for distributed arrays) that stores topology information and everything needed for communication of ghost values. This particular object is created to get partitioning information and cannot be used for real Petsc Vec and Mat objects, because they may have a different number of components.
  
  std::array<int,D> beginGlobal_;   ///< global no.s of the lower left front corner of the domain (with ghost nodes)
  std::array<node_no_t,D> localSizeWithGhosts_;     ///< local size in the coordinate directions of the local portion (including ghost nodes)
  std::array<node_no_t,D> globalSize_;    ///< global size
  std::array<int,D> nRanks_;    ///<  number of ranks in each coordinate direction that decompose the total domain
 
  std::array<std::vector<element_no_t>,D> localSizesOnRanks_;  ///< the local sizes on the ranks
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 */
template<int D, typename BasisFunctionType>
class MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>> : 
  public MeshPartitionBase
{
public:
  
  //! constructor
  MeshPartition(element_no_t globalSize, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t localSize();
  
  //! number of nodes in total
  global_no_t globalSize();
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! get an AO object
  AO &applicationOrdering();
  
  //! from a vector of global numbers remove all that are non-local
  template <typename T>
  void extractLocalNumbers(std::vector<T> &vector);

protected:
 
  element_no_t globalSize_;   ///< the global size, i.e. number of elements or nodes of the whole problem
  element_no_t localSize_;   ///< the local size, i.e. the number of elements or nodes on the local rank
};

/** partial specialization for Mesh::None, i.e. for not mesh-related partitions
 */
template<>
class MeshPartition<Mesh::None> :
  public MeshPartitionBase
{
public:
  
  //! constructor
  MeshPartition(element_no_t globalSize, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t localSize();
  
  //! number of nodes in total
  global_no_t globalSize();
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! get an AO object
  AO &applicationOrdering();
  
  //! from a vector of global numbers remove all that are non-local
  template <typename T>
  void extractLocalNumbers(std::vector<T> &vector);

protected:
 
  global_no_t globalSize_;   ///< the global size, i.e. number of elements or nodes of the whole problem
  element_no_t localSize_;   ///< the local size, i.e. the number of elements or nodes on the local rank
  global_no_t beginGlobal_;  ///< first index of the local portion
};

}  // namespace
#include "partition/01_mesh_partition_none.tpp"
#include "partition/01_mesh_partition_structured.tpp"
#include "partition/01_mesh_partition_unstructured.tpp"
