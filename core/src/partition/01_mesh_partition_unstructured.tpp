#include "partition/01_mesh_partition.h"

namespace Partition
{

// this is not yet implemented in parallel, it should however work serially
 
template<int D, typename BasisFunctionType>
MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
MeshPartition(element_no_t globalSize, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), globalSize_(globalSize)
{
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
  
  this->initializeLocalDofs();
}

template<int D, typename BasisFunctionType>
AO &MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
applicationOrdering()
{
}

//! get the local to global mapping for the current partition
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
localToGlobalMapping()
{
  PetscErrorCode ierr;
  std::vector<PetscInt> globalDofNos(localSize());
  std::iota(globalDofNos.begin(), globalDofNos.end(), 0);
  ISLocalToGlobalMapping localToGlobalMapping;
  ierr = ISLocalToGlobalMappingCreate(mpiCommunicator(), 1, localSize(), 
                                      globalDofNos.size(), PETSC_COPY_VALUES, &localToGlobalMapping); CHKERRABORT(ierr);

  return localToGlobalMapping;
}

//! number of entries in the current partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
localSize()
{
  return localSize_;
}

//! number of entries in total
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
globalSize()
{
  return globalSize_;
}
  
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
extractLocalNumbers(std::vector<T> &vector)
{
  
}
  
}  // namespace