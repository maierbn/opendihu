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

/** Composite mesh partition, this contains multiple function spaces that are connected at some edges and share nodes.
 *  Using this class it is possible to have meshes that are not
 */
template<int D, typename BasisFunctionType>
class MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>> :
  public MeshPartitionBase
{
public:

  using MeshPartitionBase::nRanks;

  //! constructor
  MeshPartition(const std::vector<std::shared_ptr<FunctionSpace<StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces,
                std::shared_ptr<RankSubset> rankSubset);

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

  const std::vector<std::shared_ptr<FunctionSpace<StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces_;

  element_no_t nElementsLocal_;   //< local number of elements of all meshes combined
  global_no_t nElementsGlobal_;   //< global number of elements of all meshes combined
  dof_no_t nDofsSharedLocal_;     //< number of non-ghost dofs that are shared, the number of distinct dofs is the total number of all subFunctionSpaces minus this value
  dof_no_t nGhostDofsSharedLocal_; //< number of ghost dofs that are shared

  std::vector<std::map<dof_no_t,std::pair<int,dof_no_t>>> sharedDofs_;   //< sharedDofs_[meshNo][dofNo] = <sameAsInMeshNo,dofNoOfThatMesh> dofs that are shared between function spaces
  std::vector<dof_no_t> nDofsLocalWithoutGhostsOnRanks_;    //< number of local dofs on every rank, this may be different than the normal value because of shared dofs

  std::vector<dof_no_t> onlyNodalDofLocalNos_;   ///< vector of local nos of the dofs, not including derivatives for Hermite
  std::vector<dof_no_t> ghostDofNosGlobalPetsc_;   ///< vector of global/petsc dof nos of the ghost dofs which are stored on the local partition
  
  dof_no_t nDofsLocalWithoutGhosts_;                       ///< number of local dofs without ghosts, cached value, the actual value is derived from nElementslocal_ and mesh type

  ISLocalToGlobalMapping localToGlobalPetscMappingDofs_;   ///< local to global mapping for dofs
};

}  // namespace

#include "partition/mesh_partition/01_mesh_partition_composite.tpp"
