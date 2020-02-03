#pragma once

#include <memory>
#include <petscdmda.h>

#include "partition/mesh_partition/00_mesh_partition_base.h"
#include "control/types.h"
#include "partition/rank_subset.h"
#include "mesh/type_traits.h"
#include "mesh/face_t.h"

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
  MeshPartition(std::array<element_no_t,MeshType::dim()> nElementsLocal, std::array<global_no_t,MeshType::dim()> nElementsGlobal,
                std::array<global_no_t,MeshType::dim()> beginElementGlobal, 
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
  element_no_t nElementsLocal(int coordinateDirection) const;
  
  //! number of elments in total
  element_no_t nElementsGlobal(int coordinateDirection) const;
  
  //! global no of first local element
  global_no_t beginElementGlobal(int coordinateDirection) const;

  //! number of nodes in total
  global_no_t nNodesGlobal() const;
  
  //! global no of first local node in the partition specified by partitionIndex or the current partition if partitionIndex == -1
  global_no_t beginNodeGlobalNatural(int coordinateDirection, int partitionIndex = -1) const;
    
  //! get the number of nodes in the global Petsc ordering that are in partitions prior to the own rank
  global_no_t beginNodeGlobalPetsc() const;

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

  //! get the global natural node no for the global coordinates of this node, this can be combined with getCoordinatesGlobal
  global_no_t getNodeNoGlobalNatural(std::array<global_no_t,MeshType::dim()> coordinatesGlobal) const;

  //! get the node no in global petsc ordering from a local node no
  global_no_t getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const;

  //! transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs
  void getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const;

  //! get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs
  global_no_t getDofNoGlobalPetsc(dof_no_t dofNoLocal) const;

  //! get the global node coordinates (x,y,z) of the node given by its local node no. This also works for ghost nodes.
  std::array<global_no_t,MeshType::dim()> getCoordinatesGlobal(node_no_t nodeNoLocal) const;

  //! get the local coordinates for a local node no, also for ghost nodes. With this method and functionSpace->getNodeNo(coordinatesLocal) it is possible to implement a global-to-local mapping.
  std::array<int,MeshType::dim()> getCoordinatesLocal(node_no_t nodeNoLocal) const;

  //! from global natural coordinates compute the local coordinates, set isOnLocalDomain to true if the node with global coordinates is in the local domain
  std::array<int,MeshType::dim()> getCoordinatesLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, bool &isOnLocalDomain) const;

  //! get the local coordinates for a local element no
  std::array<int,MeshType::dim()> getElementCoordinatesLocal(element_no_t elementNoLocal) const;

  //! get the local element no. from coordinates
  element_no_t getElementNoLocal(std::array<int,MeshType::dim()> elementCoordinates) const;

  //! get the local element no. from the global no., set isOnLocalDomain to true if the node with global coordinates is in the local domain
  element_no_t getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const;

  //! get the local node no for a global petsc node no, does not work for ghost nodes
  node_no_t getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const;

  //! get the local dof no for a global petsc dof no, does not work for ghost nodes
  dof_no_t getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const;

  //! from a vector of values of global/natural node numbers remove all that are non-local, nComponents consecutive values for each dof are assumed
  template <typename T>
  void extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents=1) const;
  
  //! from a vector of values of global/natural dofs remove all that are non-local
  template <typename T>
  void extractLocalDofsWithoutGhosts(std::vector<T> &values) const;

  //! from a vector of values of global/natural dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;

  //! get the partition index in a given coordinate direction from the rankNo
  int convertRankNoToPartitionIndex(int coordinateDirection, int rankNo);

  //! output to stream for debugging
  void output(std::ostream &stream);
  
  //! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
  //! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
  const std::vector<PetscInt> &dofNosLocal(bool onlyNodalValues=false) const;

  // use getDofNoGlobalPetsc(dofNosLocal(), ...) to get dofNosGlobalPetsc

  //! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
  void getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const;

  //! get the global dof nos of the ghost dofs in the local partition
  const std::vector<PetscInt> &ghostDofNosGlobalPetsc() const;
  
  //! Initialize the vector dofNosLocalNaturalOrdering_, this needs the functionSpace and has to be called before dofNosLocalNaturalOrdering() can be used.
  //! If the vector is already initialized by a previous call to this method, it has no effect.
  void initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>> functionSpace);

  //! Get a vector of local dof nos in local natural ordering, initializeDofNosLocalNaturalOrdering has to be called beforehand.
  const std::vector<dof_no_t> &dofNosLocalNaturalOrdering() const;

  //! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
  bool isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const;

  //! get information about neighbouring rank and boundary elements for specified face,
  //! @param neighbourRankNo: the rank of the neighbouring process that shares the face, @param nElements: Size of one-layer mesh that contains boundary elements that touch the neighbouring process
  void getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,MeshType::dim()> &nBoundaryElements, std::vector<dof_no_t> &dofNos);

  //! get the rank no of the neighbour in direction face, -1 if there is no such neighbour
  int neighbourRank(Mesh::face_t face);

  //! get the partitioning index in the coordinate direction, i.e. the no. of this rank in this direction, the total number of ranks in each direction can be retrieved by nRanks
  int ownRankPartitioningIndex(int coordinateDirection);

  //! refine the partitioning by multiplying the number of elements by refinementFactor
  void refine(std::array<int,MeshType::dim()> refinementFactor);

protected:
  
  //! initialize the values of hasFullNumberOfNodes_ variable
  void initializeHasFullNumberOfNodes();
  
  //! initialize mesh partition for a mesh with only 0 elements, 1 node and 1 dof
  void initialize1NodeMesh();

  //! initialize the value of nDofsLocalWithoutGhosts
  void setNDofsLocalWithoutGhosts();

  //! initialize the localSizesOnPartitions_ array from localSizesOnPartitions_
  void setLocalSizesOnPartitions();

  //! create the DM object for the node partitioning, such that is follows the element partitioning
  void createDmElements();
  
  //! fill the dofLocalNo vectors, onlyNodalDofLocalNos_, ghostDofNosGlobalPetsc_ and localToGlobalPetscMappingDofs_
  void createLocalDofOrderings();

  //! determine the values of ownRankPartitioningIndex_
  void setOwnRankPartitioningIndex();
  
  //! get the index in terms of partitions of the partition that contains the given node no
  std::array<int,MeshType::dim()> getPartitioningIndex(std::array<global_no_t,MeshType::dim()> nodeNoGlobalNatural) const;

  //! get the number of nodes in the global Petsc ordering that are in partitions prior to the one given by partitionIndex
  global_no_t nNodesGlobalPetscInPreviousPartitions(std::array<int,MeshType::dim()> partitionIndex) const;

  //! get the node no in global petsc ordering from global coordinates
  global_no_t getNodeNoGlobalPetsc(std::array<global_no_t,MeshType::dim()> coordinatesGlobal) const;

  std::shared_ptr<DM> dmElements_;    ///< PETSc DMDA object (data management for distributed arrays) that stores topology information and everything needed for communication of ghost values. This particular object is created to get partitioning information on the element level.
  
  std::array<global_no_t,MeshType::dim()> beginElementGlobal_;   ///< global element no.s of the lower left front corner of the domain
  std::array<node_no_t,MeshType::dim()> nElementsLocal_;     ///< local size, i.e. number of nodes in the coordinate directions of the local portion
  std::array<global_no_t,MeshType::dim()> nElementsGlobal_;    ///< global number of elements in the coodinate directions
  std::array<int,MeshType::dim()> nRanks_;    ///<  number of ranks in each coordinate direction that decompose the total domain
  std::array<int,MeshType::dim()> ownRankPartitioningIndex_;   ///< the index in terms of partitions of the own partition

  std::array<std::vector<element_no_t>,MeshType::dim()> localSizesOnPartitions_;  ///< the sizes of different partitions in each coordinate direction, i.e. localSizesOnPartitions_[0] is (width partition #0, width partition #1, ...)

  std::array<bool,MeshType::dim()> hasFullNumberOfNodes_;   ///< if the own local partition has nodes on both sides of the 1D projection at the border. This is only true at the right/top/back-most partition.
  
  std::vector<dof_no_t> onlyNodalDofLocalNos_;   ///< vector of local nos of the dofs, not including derivatives for Hermite
  std::vector<dof_no_t> ghostDofNosGlobalPetsc_;   ///< vector of global/petsc dof nos of the ghost dofs which are stored on the local partition
  
  std::vector<dof_no_t> dofNosLocalNaturalOrdering_;  ///< for every local natural number, i.e. local numbering according to coordinates, the local dof no
  dof_no_t nDofsLocalWithoutGhosts_;                       ///< number of local dofs without ghosts, cached value, the actual value is derived from nElementslocal_ and mesh type

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
  dof_no_t nDofsLocalWithGhosts() const;
  
  //! number of dofs in the local partition
  dof_no_t nDofsLocalWithoutGhosts() const;
  
  //! number of nodes in the local partition
  node_no_t nNodesLocalWithoutGhosts() const;
  
  //! return the number of local nodes for coordinateDirection == 0 and 1 otherwise (this is needed for combineFiles with paraview output writer and 3D meshes)
  node_no_t nNodesLocalWithGhosts(int coordinateDirection) const;

  //! number of nodes in the local partition
  node_no_t nNodesLocalWithGhosts() const;
  
  //! number of nodes in total
  global_no_t nNodesGlobal() const;

  //! return the number of nodes for coordinateDirection == 0 and 1 otherwise
  global_no_t nNodesGlobal(int coordinateDirection) const;

  //! number of dofs
  global_no_t nDofs() const;

  //! number of dofs in total, same as nDofs
  global_no_t nDofsGlobal() const;
  
  //! get the global element no for a local element no, this only has an effect for structured meshes, not for unstructured meshes
  global_no_t getElementNoGlobalNatural(element_no_t elementNoLocal) const;

  //! get the local to global mapping for the current partition
  ISLocalToGlobalMapping localToGlobalMappingDofs();
  
  //! from a vector of values of global node numbers remove all that are non-local, nComponents consecutive components are assumed for each dof
  template <typename T>
  void extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents=1) const;
  
  //! from a vector of values of global dofs remove all that are non-local
  void extractLocalDofsWithoutGhosts(std::vector<double> &values) const;

  //! from a vector of values of global dofs remove all that are non-local, templated type
  template <typename T>
  void extractLocalDofsWithoutGhosts(std::vector<T> &values) const;

  //! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
  void getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const;

  //! from global natural coordinates compute the local coordinates, set isOnLocalDomain to true if the node with global coordinates is in the local domain
  std::array<int,D> getCoordinatesLocal(std::array<global_no_t,D> coordinatesGlobal, bool &isOnLocalDomain) const;

  //! get the local node no for a global petsc node no, does not work for ghost nodes
  node_no_t getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const;

  //! get the local dof no for a global petsc dof no, does not work for ghost nodes
  dof_no_t getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const;

  //! get the local element no from the global element no, isOnLocalDomain is true
  element_no_t getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const;

  //! this does nothing for unstructured meshes, only for structured meshes
  void initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>> functionSpace){};

  //! output to stream for debugging
  void output(std::ostream &stream);

  //! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
  bool isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const;

  //! get the node no in global petsc ordering from a local node no
  global_no_t getNodeNoGlobalPetsc(element_no_t nodeNoLocal) const;

  //! get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs
  global_no_t getDofNoGlobalPetsc(dof_no_t dofNoLocal) const;

  //! transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs
  void getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const;

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

#include "partition/mesh_partition/01_mesh_partition_output.tpp"
#include "partition/mesh_partition/01_mesh_partition_structured.tpp"
#include "partition/mesh_partition/01_mesh_partition_structured_initialize.tpp"
#include "partition/mesh_partition/01_mesh_partition_structured_get.tpp"
#include "partition/mesh_partition/01_mesh_partition_structured_is_non_ghost.tpp"
#include "partition/mesh_partition/01_mesh_partition_structured_coordinates.tpp"
#include "partition/mesh_partition/01_mesh_partition_unstructured.tpp"
