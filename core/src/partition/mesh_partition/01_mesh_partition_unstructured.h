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

#include "partition/mesh_partition/01_mesh_partition_unstructured.tpp"
