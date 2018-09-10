#include "partition/01_mesh_partition.h"

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
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesLocalWithGhosts() const
{
  return nNodes_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nDofsLocalWithGhosts() const
{
  return nDofs_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nDofsLocalWithoutGhosts() const
{
  return nDofs_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
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
nDofs() const
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
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
output(std::ostream &stream)
{
  stream << "MeshPartition<Unstructured>, nElements_: " << nElements_ << ", nNodes_: " << nNodes_ << ", nDofs_: " << nDofs_;
}

}  // namespace
