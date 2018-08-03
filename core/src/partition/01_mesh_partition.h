#pragma once

#include <memory>
#include <petscdmda.h>

#include "partition/00_mesh_partition_base.h"
#include "control/types.h"
#include "partition/rank_subset.h"
#include "mesh/type_traits.h"

// forward declaration
namespace BasisOnMesh 
{
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh; 
}

namespace Partition
{
 
/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: ghost elements
 */
template<typename BasisOnMeshType, typename DummyForTraits = typename BasisOnMeshType::Mesh>
class MeshPartition
{
};

/** partial specialization for structured meshes */
template<typename MeshType, typename BasisFunctionType>
class MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>> :
  public MeshPartitionBase
{
public:
 
  //! constructor, determine the decomposition by PETSc
  MeshPartition(std::array<global_no_t,MeshType::dim()> globalSize, std::shared_ptr<RankSubset> rankSubset);
 
  //! constructor from prescribed partition
  MeshPartition(std::array<node_no_t,MeshType::dim()> localSize, std::array<global_no_t,MeshType::dim()> globalSize,
                std::array<int,MeshType::dim()> beginGlobal, 
                std::array<int,MeshType::dim()> nRanks, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of ranks in a coordinate direction
  int nRanks(int coordinateDirection);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t localSize();
  
  //! number of entries in the given coordinate direction in the current partition
  element_no_t localSize(int coordinateDirection);

  //! global no of first local number in current partition
  global_no_t beginGlobal(int coordinateDirection);
  
  //! global no of first local node in current partition
  global_no_t beginNodeGlobal(int coordinateDirection);
  
  //! global no of one after last local number in current partition
  global_no_t endGlobal(int coordinateDirection);
  
  //! number of nodes in total
  global_no_t globalSize();
  
  //! number of nodes in total, in the given coordinate direction 
  global_no_t globalSize(int coordinateDirection);
  
  //! get if there are nodes on both borders in the given coordinate direction
  //! this is the case if the local partition touches the right/top/back border
  bool hasFullNumberOfNodes(int coordinateDirection);
  
  //! get a vector with the local sizes on every rank
  std::vector<element_no_t> &localSizesOnRanks(int coordinateDirection);
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! from a vector of values of global node numbers remove all that are non-local
  template <typename T>
  void extractLocalNodes(std::vector<T> &vector);
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofs(std::vector<double> &values);
  
protected:
 
  DM dm_;    ///< PETSc DMDA object (data management for distributed arrays) that stores topology information and everything needed for communication of ghost values. This particular object is created to get partitioning information and cannot be used for real Petsc Vec and Mat objects, because they may have a different number of components.
  
  std::array<int,MeshType::dim()> beginGlobal_;   ///< global no.s of the lower left front corner of the domain (with ghost nodes)
  std::array<node_no_t,MeshType::dim()> localSizeWithGhosts_;     ///< local size in the coordinate directions of the local portion (including ghost nodes)
  std::array<global_no_t,MeshType::dim()> globalSize_;    ///< global size
  std::array<int,MeshType::dim()> nRanks_;    ///<  number of ranks in each coordinate direction that decompose the total domain
 
  std::array<std::vector<element_no_t>,MeshType::dim()> localSizesOnRanks_;  ///< the local sizes on the ranks
  std::array<bool,MeshType::dim()> hasFullNumberOfNodes_;   ///< if the own local partition has nodes on both sides of the 1D projection at the border. This is only true at the right/top/back-most partition.
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
  MeshPartition(global_no_t globalSize, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t localSize();
  
  //! number of nodes in total
  global_no_t globalSize();
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! from a vector of values of global node numbers remove all that are non-local
  template <typename T>
  void extractLocalNodes(std::vector<T> &vector);
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofs(std::vector<double> &values);
  
protected:
 
  global_no_t globalSize_;   ///< the global size, i.e. number of elements or nodes of the whole problem
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
  MeshPartition(global_no_t globalSize, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t localSize();
  
  //! number of nodes in total
  global_no_t globalSize();
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! from a vector of values of global node numbers remove all that are non-local
  template <typename T>
  void extractLocalNodes(std::vector<T> &vector);
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofs(std::vector<double> &values);
  
  
protected:
 
  global_no_t globalSize_;   ///< the global size, i.e. number of elements or nodes of the whole problem
  element_no_t localSize_;   ///< the local size, i.e. the number of elements or nodes on the local rank
  global_no_t beginGlobal_;  ///< first index of the local portion
};

}  // namespace
#include "partition/01_mesh_partition_none.tpp"
#include "partition/01_mesh_partition_structured.tpp"
#include "partition/01_mesh_partition_unstructured.tpp"
