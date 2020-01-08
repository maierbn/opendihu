#include "partition/mesh_partition/01_mesh_partition.h"

namespace Partition
{

// this is not yet implemented in parallel, it should however work serially
 
template<int D, typename BasisFunctionType>
MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
MeshPartition(global_no_t nElementsGlobal, global_no_t nNodesGlobal, global_no_t nDofsGlobal, std::shared_ptr<RankSubset> rankSubset):
  MeshPartitionBase(rankSubset), nElements_(nElementsGlobal), nNodes_(nNodesGlobal), nDofs_(nDofsGlobal)
{
  global_no_t nDofsLocal = nDofsGlobal;
  this->createLocalDofOrderings(nDofsLocal);
}

//! get the local to global mapping for the current partition
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
localToGlobalMappingDofs()
{
  PetscErrorCode ierr;
  std::vector<PetscInt> globalDofNos(nDofs_);
  std::iota(globalDofNos.begin(), globalDofNos.end(), 0);
  ISLocalToGlobalMapping localToGlobalMapping;
  ierr = ISLocalToGlobalMappingCreate(mpiCommunicator(), 1, nDofs_, 
                                      globalDofNos.data(), PETSC_COPY_VALUES, &localToGlobalMapping); CHKERRABORT(mpiCommunicator(),ierr);

  return localToGlobalMapping;
}

//! number of entries in the current partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nElementsLocal() const
{
  return nElements_;
}

//! number of entries in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nElementsGlobal() const
{
  return nElements_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesLocalWithGhosts() const
{
  return nNodes_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nDofsLocalWithGhosts() const
{
  return nDofs_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nDofsLocalWithoutGhosts() const
{
  return nDofs_;
}


//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesLocalWithGhosts(int coordinateDirection) const
{
  if (coordinateDirection == 0)
    return nNodes_;
  return 1;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesLocalWithoutGhosts() const
{
  return nNodes_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesGlobal() const
{
  return nNodes_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesGlobal(int coordinateDirection) const
{
  if (coordinateDirection == 0)
    return nNodes_;
  return 1;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nDofs() const
{
  return nDofs_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nDofsGlobal() const
{
  return nDofs_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getElementNoGlobalNatural(element_no_t elementNoLocal) const
{
  return (global_no_t)(elementNoLocal);
}
  
//! get the local element no. from the global no., set isOnLocalDomain to true if the node with global coordinates is in the local domain
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const
{
  isOnLocalDomain = true;
  return elementNoGlobalPetsc;
}

template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents) const
{
  
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<double> &vector) const
{

}


template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<T> &vector) const
{

}

template<int D, typename BasisFunctionType>
std::array<int,D> MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getCoordinatesLocal(std::array<global_no_t,D> coordinatesGlobal, bool &isOnLocalDomain) const
{
  // because we have no parallel execution, global coordinates are the same as local coordinates and isOnLocalDomain is always true
  isOnLocalDomain = true;
  std::array<int,D> coordinatesLocal(coordinatesGlobal.begin(), coordinatesGlobal.end());
  return coordinatesLocal;
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const
{
  dofNosGlobalNatural.resize(nDofsLocalWithoutGhosts());
  std::iota(dofNosGlobalNatural.begin(), dofNosGlobalNatural.end(), 0);
}

template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const
{
  isLocal = true;
  return (node_no_t)nodeNoGlobalPetsc;
}

template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const
{
  isLocal = true;
  return (dof_no_t)dofNoGlobalPetsc;
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
output(std::ostream &stream)
{
  stream << "MeshPartition<Unstructured>, nElements_: " << nElements_ << ", nNodes_: " << nNodes_ << ", nDofs_: " << nDofs_;
}

//! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
template<int D, typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const
{
  return true;
}

//! get the node no in global petsc ordering from a local node no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const
{
  return (global_no_t)nodeNoLocal;
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const
{
  dofNosGlobalPetsc.assign(dofNosLocal.begin(), dofNosLocal.end());
}

template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getDofNoGlobalPetsc(dof_no_t dofNoLocal) const
{
  return (global_no_t)dofNoLocal;
}

}  // namespace
