#pragma once

#include <memory>
#include <petscdmda.h>

#include "utility/vector_operators.h"
#include "partition/mesh_partition/00_mesh_partition_base.h"
#include "partition/mesh_partition/01_mesh_partition.h"
#include "control/types.h"
#include "partition/rank_subset.h"
#include "mesh/type_traits.h"
#include "mesh/face_t.h"
#include "mesh/structured_deformable.h"

// forward declaration
namespace FunctionSpace 
{
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace; 
}

namespace Partition
{

/** Composite mesh partition, this contains multiple function spaces with mesh type Mesh::StructuredDeformableOfDimension<D>
 *  that are connected at some edges and share nodes.
 *  Using this class it is possible to have meshes that are not cuboid but multiple cuboids connected together.
 *  There are shared nodes that will only appear once in the numbering and therefore share the dofs.
 *  In a parallel setting, shared nodes have to be on the same subdomain.
 *
 *  The composite mesh appears like a normal structured mesh except it does not have the coordinate dependent methods.
 *  Internally, a composite numbering is created that avoids duplicate=shared nodes. For every duplicate node, the dofs are only included once.
 *  The initialization of this numbering involves communication of ghost node no with the neighbours.
 */
template<int D, typename BasisFunctionType>
class MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>> :
  public MeshPartitionBase
{
public:

  using MeshPartitionBase::nRanks;

  typedef FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! constructor
  MeshPartition(const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces,
                std::shared_ptr<RankSubset> rankSubset);

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
  
  //! number of nodes in total
  global_no_t nNodesGlobal() const;
  
  //! get the number of nodes in the global Petsc ordering that are in partitions prior to the own rank
  global_no_t beginNodeGlobalPetsc() const;

  //! get the local to global mapping for the current partition, for the dof numbering
  ISLocalToGlobalMapping localToGlobalMappingDofs();
  
  //! get the global natural element no for a local element no, this is used for parsing elemental boundary conditions
  //! the global natural ordering is here defined to be the global natural ordering of all submeshes concatenated
  global_no_t getElementNoGlobalNatural(element_no_t elementNoLocal) const;

  //! get the node no in global petsc ordering from a local node no (works also for ghosts)
  global_no_t getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const;

  //! get the node no in a composite global natural ordering where the natural orders of the submeshes are concatenated, call this method from the function space to be compatible with structured meshes!
  global_no_t getNodeNoGlobalNatural(element_no_t elementNoLocal, int nodeIndex) const;

  //! transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs
  void getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const;

  //! get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs
  global_no_t getDofNoGlobalPetsc(dof_no_t dofNoLocal) const;

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

  //! output to stream for debugging
  void output(std::ostream &stream);
  
  //! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
  //! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
  const std::vector<PetscInt> &dofNosLocal(bool onlyNodalValues=false) const;

  // use getDofNoGlobalPetsc(dofNosLocal(), ...) to get dofNosGlobalPetsc

  //! get the global dof nos of the ghost dofs in the local partition
  const std::vector<PetscInt> &ghostDofNosGlobalPetsc() const;
  
  //! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
  bool isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const;

  //! get information about neighbouring rank and boundary elements for specified face,
  //! @param neighbourRankNo: the rank of the neighbouring process that shares the face, @param nElements: Size of one-layer mesh that contains boundary elements that touch the neighbouring process
  void getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,D> &nBoundaryElements, std::vector<dof_no_t> &dofNos);

  //! get the rank no of the neighbour in direction face, -1 if there is no such neighbour
  int neighbourRank(Mesh::face_t face);

  //! from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
  void getSubMeshNoAndElementNoLocal(element_no_t elementNoLocal, int &subMeshNo, element_no_t &elementOnMeshNoLocal) const;

  //! from a local node no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
  void getSubMeshNoAndNodeNoLocal(node_no_t nodeNoLocal, int &subMeshNo, node_no_t &nodeOnMeshNoLocal) const;

  //! for the local node no in the composite numbering return all sub meshes and the corresponding local node nos in non-composite numbering of this node. This may be multiple if the node is shared.o
  void getSubMeshesWithNodes(node_no_t nodeNoLocal, std::vector<std::pair<int,node_no_t>> &subMeshesWithNodes) const;

  //! from the submesh no and the local node no in the submesh numbering get the local node no in the composite numbering
  node_no_t getNodeNoLocalFromSubmesh(int subMeshNo, node_no_t nodeNoDuplicateOnSubmesh) const;

  //! get a string with all information, this is used in the regression tests (unit test) to compare it to a reference string
  std::string getString();

protected:

  //! initialize removedSharedNodes_, find nodes that are the same on multiple meshes
  void initializeSharedNodes();

  //! initialize the numberings
  void initializeGhostNodeNos();

  //! fill the dofLocalNo vectors, onlyNodalDofLocalNos_, ghostDofNosGlobalPetsc_ and localToGlobalPetscMappingDofs_
  void createLocalDofOrderings();

  int nSubMeshes_;                //< number of sub function spaces, = subFunctionSpaces_.size()
  const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces_;

  std::vector<std::map<node_no_t,std::pair<int,node_no_t>>> removedSharedNodes_;   //< removedSharedNodes_[meshNo][nodeNo] = <sameAsInMeshNo,nodeNoOfThatMesh> nodes that are shared between function spaces, they appear only once in the second function space and are removed there (not included in the composite mapping)

  element_no_t nElementsLocal_;   //< local number of elements of all meshes combined
  global_no_t nElementsGlobal_;   //< global number of elements of all meshes combined
  global_no_t elementNoGlobalBegin_; //< first global element no of the local domain

  node_no_t nNodesSharedLocal_;      //< number of non-ghost nodes that are shared/removed, the number of distinct nodes is the total number of all subFunctionSpaces minus this value
  node_no_t nGhostNodesSharedLocal_; //< number of ghost nodes that are shared/removed

  // -------------- everything below is initialized by initializeGhostNodeNos()
  std::vector<int> nRemovedNodesNonGhost_;            //< for every mesh the number of duplicate (non-ghost) nodes that are thus not included in the mesh
  std::vector<int> nNonDuplicateNodesWithoutGhosts_;  //< the local without ghosts number of nodes in each submesh
  node_no_t nNodesLocalWithoutGhosts_;                //< number of local nodes without ghosts in the total mesh in the duplicate-free numbering
  node_no_t nNodesLocalWithGhosts_;                   //< number of local nodes with ghosts in the total mesh in the duplicate-free numbering
  global_no_t nNodesGlobal_;                          //< the global number of nodes on all submeshes
  global_no_t nonDuplicateNodeNoGlobalBegin_;         //< the first no in the duplicate-free global numbering

  struct NodesRequest
  {
    std::vector<global_no_t> nodeNosGlobalPetsc;      //< global node no
    std::vector<node_no_t> nodeNosLocal;              //< local node no on own rank
  };

  std::map<int, std::vector<NodesRequest>> requestNodesFromRanks_;    //< requestNodesFromRanks_[rankNo][subMeshNo].nodeNosGlobalPetsc, for some other ranks which nodes are requested from them, for each submesh

  // mappings from local numbering in every submesh to the composite numbering
  std::vector<std::vector<node_no_t>> meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_;   //< mapping from submesh no and local node no to the composite numbering used for the whole mesh, for local nodes with ghosts, -1 for removed nodes
  std::vector<std::vector<node_no_t>> meshAndNodeNoLocalToNodeNoNonDuplicateLocal_;    //< mapping from submesh no and local node no to the local number of the composite node numbering, also for ghost nodes
  std::vector<std::vector<bool>> isDuplicate_;                //< for every local node no, if the node has a prescribed Dirichlet BC value
  std::vector<std::pair<int,node_no_t>> nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_;   //< mapping from non-duplicate local number to submesh no and local node no on the submesh

  std::vector<PetscInt> nonDuplicateGhostNodeNosGlobal_;  //< duplicate-free ghost nodes in duplicate-free global indexing, needed to create Petsc Vecs

  // ------------- everything below is initialized in createLocalDofOrderings()
  std::vector<node_no_t> onlyNodalDofLocalNos_;       //< vector of local dofs of the nodes, not including derivatives for Hermite
  std::vector<PetscInt> ghostDofNosGlobalPetsc_;      //< vector of global/petsc dof nos of the ghost nodes which are stored on the local partition

  ISLocalToGlobalMapping localToGlobalPetscMappingDofs_;   //< local to global mapping for nodes

};

}  // namespace

#include "partition/mesh_partition/01_mesh_partition_composite.tpp"
#include "partition/mesh_partition/01_mesh_partition_composite_initialize.tpp"
