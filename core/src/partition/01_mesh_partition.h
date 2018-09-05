#pragma once

#include <memory>
#include <petscdmda.h>

#include "partition/00_mesh_partition_base.h"
#include "control/types.h"
#include "partition/rank_subset.h"
#include "mesh/type_traits.h"

// forward declaration
namespace FunctionSpace 
{
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace; 
}

namespace Partition
{
 
/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: starting with 0, up to total number including ghost elements
 */
template<typename FunctionSpaceType, typename DummyForTraits = typename FunctionSpaceType::Mesh>
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
 *  for (node_no_t localNodeNo = 0; localNodeNo < meshPartition_->nNodesLocalWithoutGhosts(); localNodeNo++)
 *  {
 *    for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++)
 *    {
 *      dof_no_t dofLocalNo = localNodeNo*nDofsPerNode + dofIndex;
 *    }
 *  }
 */
template<typename MeshType, typename BasisFunctionType>
class MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>> :
  public MeshPartitionBase
{
public:

  using MeshPartitionBase::nRanks;
 
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
  
  //! number of dofs in total
  global_no_t nDofsGlobal() const;

  //! number of nodes in the local partition
  node_no_t nNodesLocalWithGhosts() const;
  
  //! number of nodes in the local partition
  node_no_t nNodesLocalWithoutGhosts() const;
  
  //! number of nodes in the local partition specified by partitionIndex or the current partition if partitionIndex == -1
  node_no_t nNodesLocalWithGhosts(int coordinateDirection, int partitionIndex = -1) const;
  
  //! number of nodes in the partition specified by partitionIndex or the current partition if partitionIndex == -1
  node_no_t nNodesLocalWithoutGhosts(int coordinateDirection, int partitionIndex = -1) const;
  
  //! number of elments in the local partition
  node_no_t nElementsLocal(int coordinateDirection) const;
  
  //! number of elments in total
  node_no_t nElementsGlobal(int coordinateDirection) const;
  
  //! number of nodes in total
  global_no_t nNodesGlobal() const;
  
  //! global no of first local node in the partition specified by partitionIndex or the current partition if partitionIndex == -1
  global_no_t beginNodeGlobalNatural(int coordinateDirection, int partitionIndex = -1) const;
    
  //! number of nodes in total
  global_no_t nNodesGlobal(int coordinateDirection) const;
  
  //! get if there are nodes on both borders in the given coordinate direction
  //! this is the case if the partition touches the right/top/back border
  //! Consider the partition specified by partitionIndex or the current partition if partitionIndex == -1.
  bool hasFullNumberOfNodes(int coordinateDirection, int partitionIndex = -1) const;
  
  //! get a vector with the local sizes on every rank
  const std::vector<element_no_t> &localSizesOnRanks(int coordinateDirection) const;
  
  //! get the local to global mapping for the current partition, for the dof numbering
  ISLocalToGlobalMapping localToGlobalMappingDofs();
  
  //! get the global natural element no for a local element no
  global_no_t getElementNoGlobalNatural(element_no_t elementNoLocal) const;

  //! get the global natural node no for the global coordinates of this node, this can be combined with getNodeNoGlobalCoordinates
  global_no_t getNodeNoGlobalNatural(std::array<int,MeshType::dim()> coordinates) const;

  //! get the global node coordinates (x,y,z) of the node given by its local node no. This also works for ghost nodes.
  std::array<int,MeshType::dim()> getNodeNoGlobalCoordinates(node_no_t nodeNoLocal) const;

  //! from a vector of values of global/natural node numbers remove all that are non-local, nComponents consecutive values for each dof are assumed
  template <typename T>
  void extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents=1) const;
  
  //! from a vector of values of global/natural dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;
  
  //! get the partition index in a given coordinate direction from the rankNo
  int convertRankNoToPartitionIndex(int coordinateDirection, int rankNo);

  //! output to stream for debugging
  void output(std::ostream &stream);
  
  //! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
  //! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
  const std::vector<PetscInt> &dofNosLocal(bool onlyNodalValues=false) const;
  
  //! get the global dof nos of the ghost dofs in the local partition
  const std::vector<PetscInt> &ghostDofNosGlobalPetsc() const;
  
  //! Initialize the vector dofNosLocalNaturalOrdering_, this needs the functionSpace and has to be called before dofNosLocalNaturalOrdering() can be used.
  //! If the vector is already initialized by a previous call to this method, it has no effect.
  void initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>> functionSpace);

  //! Get a vector of local dof nos in local natural ordering, initializeDofNosLocalNaturalOrdering has to be called beforehand.
  const std::vector<dof_no_t> &dofNosLocalNaturalOrdering() const;

protected:
  
  //! initialize the values of hasFullNumberOfNodes_ variable
  void initializeHasFullNumberOfNodes();
  
  //! create the DM object for the node partitioning, such that is follows the element partitioning
  void createDmElements();
  
  //! fill the dofLocalNo vectors
  void createLocalDofOrderings();

  //! determine the values of ownRankPartitioningIndex_
  void setOwnRankPartitioningIndex();
  
  //! get the index in terms of partitions of the partition that contains the given node no
  std::array<int,MeshType::dim()> getPartitioningIndex(std::array<global_no_t,MeshType::dim()> nodeNoGlobalNatural);

  //! get the number of nodes in the globalPetsc ordering that in partitions prior to the one given by partitionIndex
  global_no_t nNodesGlobalPetscInPreviousPartitions(std::array<int,MeshType::dim()> partitionIndex) const;

  std::shared_ptr<DM> dmElements_;    ///< PETSc DMDA object (data management for distributed arrays) that stores topology information and everything needed for communication of ghost values. This particular object is created to get partitioning information on the element level.
  
  std::array<int,MeshType::dim()> beginElementGlobal_;   ///< global element no.s of the lower left front corner of the domain
  std::array<node_no_t,MeshType::dim()> nElementsLocal_;     ///< local size, i.e. number of nodes in the coordinate directions of the local portion (including ghost elements)
  std::array<global_no_t,MeshType::dim()> nElementsGlobal_;    ///< global number of elements in the coodinate directions
  std::array<int,MeshType::dim()> nRanks_;    ///<  number of ranks in each coordinate direction that decompose the total domain
  std::array<int,MeshType::dim()> ownRankPartitioningIndex_;   ///< the index in terms of partitions of the own partition

  std::array<std::vector<element_no_t>,MeshType::dim()> localSizesOnRanks_;  ///< the sizes of different partitions in each coordinate direction, i.e. localSizesOnRanks_[0] is (width partition #0, width partition #1, ...)

  std::array<bool,MeshType::dim()> hasFullNumberOfNodes_;   ///< if the own local partition has nodes on both sides of the 1D projection at the border. This is only true at the right/top/back-most partition.
  
  std::vector<dof_no_t> onlyNodalDofLocalNos_;   ///< vector of local nos of the dofs, not including derivatives for Hermite
  std::vector<dof_no_t> ghostDofNosGlobalPetsc_;   ///< vector of global/petsc dof nos of the ghost dofs which are stored on the local partition
  
  std::vector<dof_no_t> dofNosLocalNaturalOrdering_;  ///< for every local natural number, i.e. local numbering according to coordinates, the local dof no

  ISLocalToGlobalMapping localToGlobalPetscMappingDofs_;   ///< local to global mapping for dofs
};

/** Partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 *  Unstructured meshes should work for non-parallel execution, as usual.
 */
template<int D, typename BasisFunctionType>
class MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>> : 
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
  global_no_t getElementNoGlobalNatural(element_no_t elementNoLocal) const;

  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMappingDofs();
  
  //! from a vector of values of global node numbers remove all that are non-local, nComponents consecutive components are assumed for each dof
  template <typename T>
  void extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents=1) const;
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;

  //! this does nothing for unstructured meshes, only for structured meshes
  void initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>> functionSpace){};

  //! output to stream for debugging
  void output(std::ostream &stream);
  
protected:
 
  global_no_t nElements_;   ///< the global size, i.e. number of elements of the whole problem
  global_no_t nNodes_;   ///< the global size, i.e. the number of nodes of the whole problem
  global_no_t nDofs_;    ///< the number of dofs
};

}  // namespace

//! output mesh partition
template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, Partition::MeshPartition<FunctionSpaceType> meshPartition);
template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition);

//! output local to global mapping
std::ostream &operator<<(std::ostream &stream, ISLocalToGlobalMapping localToGlobalMapping);
std::ostream &operator<<(std::ostream &stream, std::shared_ptr<ISLocalToGlobalMapping> localToGlobalMapping);

#include "partition/01_mesh_partition_output.tpp"
#include "partition/01_mesh_partition_structured.tpp"
#include "partition/01_mesh_partition_unstructured.tpp"
