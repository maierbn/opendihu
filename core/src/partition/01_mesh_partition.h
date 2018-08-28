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
 *  Local numbering: starting with 0, up to total number including ghost elements
 */
template<typename BasisOnMeshType, typename DummyForTraits = typename BasisOnMeshType::Mesh>
class MeshPartition
{
};

/** Partial specialization for structured meshes.
 *  The items under consideration are elements, i.e. localSize, etc. refer to number of elements.
 * 
 *  You can iterate over non-ghost and ghost nodes together using
 *  for (node_no_t localNodeNo = 0; localNodeNo < meshPartition_->nNodesLocalWithGhosts(); localNodeNo++)
 *  {
 *    for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++)
 *    {
 *      dof_no_t dofLocalNo = localNodeNo*nDofsPerNode + dofIndex;
 *    }
 *  }
 * 
 *  To iterate only over non-ghost dofs:
 *  for (std::vector<dof_no_t>::const_iterator localDof = meshPartition_->nonGhostDofsBegin(); localDof != meshPartition_->nonGhostDofsEnd(); localDof++)
 *  {
 *     
 *  }
 */
template<typename MeshType, typename BasisFunctionType>
class MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>> :
  public MeshPartitionBase
{
public:
 
  //! constructor, determine the decomposition by PETSc
  MeshPartition(std::array<global_no_t,MeshType::dim()> nElementsGlobal, std::shared_ptr<RankSubset> rankSubset);
 
  //! constructor from prescribed partition
  MeshPartition(std::array<node_no_t,MeshType::dim()> nElementsLocal, std::array<global_no_t,MeshType::dim()> nElementsGlobal,
                std::array<int,MeshType::dim()> beginElementGlobal, 
                std::array<int,MeshType::dim()> nRanks, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of ranks in a coordinate direction
  int nRanks(int coordinateDirection) const;
  
  //! number of elements in the current partition
  element_no_t nElementsLocal() const;
  
  //! number of elements in total
  global_no_t nElementsGlobal() const;
  
  //! number of dofs in the local partition
  dof_no_t nDofsLocalWithGhosts() const;
  
  //! number of dofs in the local partition, without ghosts
  dof_no_t nDofsLocalWithoutGhosts() const;
  
  //! number of nodes in the local partition
  node_no_t nNodesLocalWithGhosts() const;
  
  //! number of nodes in the local partition
  node_no_t nNodesLocalWithoutGhosts() const;
  
  //! number of nodes in the local partition
  node_no_t nNodesLocalWithGhosts(int coordinateDirection) const;
  
  //! number of nodes in the local partition
  node_no_t nNodesLocalWithoutGhosts(int coordinateDirection) const;
  
  //! number of elments in the local partition
  node_no_t nElementsLocal(int coordinateDirection) const;
  
  //! number of elments in total
  node_no_t nElementsGlobal(int coordinateDirection) const;
  
  //! number of nodes in total
  global_no_t nNodesGlobal() const;
  
  //! global no of first local node in current partition
  node_no_t beginNodeGlobal(int coordinateDirection) const;
    
  //! number of nodes in total
  global_no_t nNodesGlobal(int coordinateDirection) const;
  
  //! get if there are nodes on both borders in the given coordinate direction
  //! this is the case if the local partition touches the right/top/back border
  bool hasFullNumberOfNodes(int coordinateDirection) const;
  
  //! get a vector with the local sizes on every rank
  const std::vector<element_no_t> &localSizesOnRanks(int coordinateDirection) const;
  
  //! get the local to global mapping for the current partition, for the dof numbering
  ISLocalToGlobalMapping localToGlobalMappingDofs();
  
  //! get the global element no for a local element no
  global_no_t getElementNoGlobal(element_no_t elementNoLocal) const;

  //! from a vector of values of global node numbers remove all that are non-local
  template <typename T>
  void extractLocalNodes(std::vector<T> &vector) const;
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;
  
  //! output to stream for debugging
  void output(std::ostream &stream);
  
  //! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
  //! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
  const std::vector<PetscInt> &dofNosLocal(bool onlyNodalValues=false) const;
  
  //! get the global dof nos of the ghost dofs in the local partition
  const std::vector<PetscInt> &ghostDofGlobalNos() const;
  
protected:
  
  //! initialize the values of hasFullNumberOfNodes_ variable
  void initializeHasFullNumberOfNodes();
  
  //! create the DM object for the node partitioning, such that is follows the element partitioning
  void createDmElements();
  
  //! fill the dofLocalNo vectors
  void createLocalDofOrderings();
  
  std::shared_ptr<DM> dmElements_;    ///< PETSc DMDA object (data management for distributed arrays) that stores topology information and everything needed for communication of ghost values. This particular object is created to get partitioning information on the element level.
  
  std::array<int,MeshType::dim()> beginElementGlobal_;   ///< global element no.s of the lower left front corner of the domain
  std::array<node_no_t,MeshType::dim()> nElementsLocal_;     ///< local size, i.e. number of nodes in the coordinate directions of the local portion (including ghost elements)
  std::array<global_no_t,MeshType::dim()> nElementsGlobal_;    ///< global number of elements in the coodinate directions
  std::array<int,MeshType::dim()> nRanks_;    ///<  number of ranks in each coordinate direction that decompose the total domain
 
  std::array<std::vector<element_no_t>,MeshType::dim()> localSizesOnRanks_;  ///< the local sizes on the ranks
  std::array<bool,MeshType::dim()> hasFullNumberOfNodes_;   ///< if the own local partition has nodes on both sides of the 1D projection at the border. This is only true at the right/top/back-most partition.
  
  std::vector<dof_no_t> onlyNodalDofLocalNos_;   ///< vector of local nos of the dofs, not including derivatives for Hermite
  std::vector<dof_no_t> ghostDofGlobalNos_;   ///< vector of global dof nos of the ghost dofs which are stored on the local partition
  
  ISLocalToGlobalMapping localToGlobalMapping_;   ///< local to global mapping for dofs
};

/** Partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 *  Unstructured meshes should work for non-parallel execution, as usual.
 */
template<int D, typename BasisFunctionType>
class MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>> : 
  public MeshPartitionBase
{
public:
  
  //! constructor
  MeshPartition(global_no_t nElementsGlobal, global_no_t nNodesGlobal, global_no_t nDofsGlobal, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of elements in the local partition
  element_no_t nElementsLocal() const;
  
  //! number of elements in total
  global_no_t nElementsGlobal() const;
  
  //! number of dofs in the local partition
  element_no_t nDofsLocalWithGhosts() const;
  
  //! number of dofs in the local partition
  element_no_t nDofsLocalWithoutGhosts() const;
  
  //! number of nodes in the local partition
  element_no_t nNodesLocalWithoutGhosts() const;
  
  //! number of nodes in the local partition
  element_no_t nNodesLocalWithGhosts() const;
  
  //! number of nodes in total
  global_no_t nNodesGlobal() const;
  
  //! number of dofs
  global_no_t nDofs() const;
  
  //! get the global element no for a local element no, this only has an effect for structured meshes, not for unstructured meshes
  global_no_t getElementNoGlobal(element_no_t elementNoLocal) const;

  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMappingDofs();
  
  //! from a vector of values of global node numbers remove all that are non-local
  template <typename T>
  void extractLocalNodes(std::vector<T> &vector) const;
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;
  
  //! output to stream for debugging
  void output(std::ostream &stream);
  
protected:
 
  global_no_t nElements_;   ///< the global size, i.e. number of elements of the whole problem
  global_no_t nNodes_;   ///< the global size, i.e. the number of nodes of the whole problem
  global_no_t nDofs_;    ///< the number of dofs
};

/** Partial specialization for Mesh::None, i.e. for not mesh-related partitions.
 *  The items of this partition are not necessarily elements or nodes.
 * TODO: delete
 */
template<>
class MeshPartition<Mesh::None> :
  public MeshPartitionBase
{
public:
  
  //! constructor
  MeshPartition(global_no_t globalSize, std::shared_ptr<RankSubset> rankSubset);
  
  //! number of entries in the current partition (this usually refers to the elements)
  element_no_t nElementsLocal() const;
  
  //! number of nodes in total
  global_no_t nElementsGlobal() const;
  
  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMapping();
  
  //! from a vector of values of global node numbers remove all that are non-local
  template <typename T>
  void extractLocalNodes(std::vector<T> &vector) const;
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;
  
  //! output to stream for debugging
  void output(std::ostream &stream);
  
protected:
 
  global_no_t globalSize_;   ///< the global size, i.e. number of elements or nodes of the whole problem
  element_no_t localSize_;   ///< the local size, i.e. the number of elements or nodes on the local rank
  global_no_t beginGlobal_;  ///< first index of the local portion
};

}  // namespace

//! output mesh partition
template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, Partition::MeshPartition<BasisOnMeshType> meshPartition);
template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition);

//! output local to global mapping
std::ostream &operator<<(std::ostream &stream, ISLocalToGlobalMapping localToGlobalMapping);
std::ostream &operator<<(std::ostream &stream, std::shared_ptr<ISLocalToGlobalMapping> localToGlobalMapping);

#include "partition/01_mesh_partition_none.tpp"
#include "partition/01_mesh_partition_structured.tpp"
#include "partition/01_mesh_partition_unstructured.tpp"
