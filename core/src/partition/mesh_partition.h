#pragma once

#include <memory>
#include <petscdmda.h>

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
class MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>
{
public:
 
  //! constructor, determine the decomposition by PETSc
  MeshPartition(std::array<node_no_t,D> globalSize, std::shared_ptr<RankSubset> rankSubset);
 
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
  
  //! get a vector with the local sizes on every rank
  std::vector<element_no_t> &localSizesOnRanks(int coordinateDirection);
  
  //! get the MPI communicator that is needed for the work portion
  MPI_Comm mpiCommunicator();
  
  //! get an AO object
  AO &applicationOrdering();
  
  //! from a vector of global numbers remove all that are non-local
  template <typename T>
  void extractLocalNumbers(std::vector<T> &vector);
  
protected:
 
  DM dm_;    ///< PETSc DMDA object (data management for distributed arrays) that stores topology information and everything needed for communication of ghost values. This particular object is created to get partitioning information and cannot be used for real Petsc Vec and Mat objects, because they may have a different number of components.
  
  std::array<global_no_t,D> beginGlobal_;   ///< global no.s of the lower left front corner of the domain (with ghost nodes)
  std::array<node_no_t,D> localSizeWithGhosts_;     ///< local size in the coordinate directions of the local portion (including ghost nodes)
  std::array<node_no_t,D> globalSize_;    ///< global size
  std::array<int,D> nRanks_;    ///<  number of ranks in each coordinate direction that decompose the total domain
 
  std::shared_ptr<RankSubset> rankSubset_;  ///< the set of ranks that compute something where this partition is a part of, also holds the MPI communciator
  std::array<std::vector<element_no_t>,D> localSizesOnRanks_;  ///< the local sizes on the ranks
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 */
template<int D, typename BasisFunctionType>
class MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>
{
public:
  
};

}  // namespace
#include "partition/mesh_partition_structured.tpp"
