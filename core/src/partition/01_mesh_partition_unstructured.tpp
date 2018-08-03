#include "partition/01_mesh_partition.h"

namespace Partition
{

// this is not yet implemented in parallel, it should however work serially
 
template<int D, typename BasisFunctionType>
MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
MeshPartition(global_no_t globalSize, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), globalSize_(globalSize)
{
  this->initializeLocalDofs();
}

//! get the local to global mapping for the current partition
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
localToGlobalMapping()
{
  PetscErrorCode ierr;
  std::vector<PetscInt> globalDofNos(localSize());
  std::iota(globalDofNos.begin(), globalDofNos.end(), 0);
  ISLocalToGlobalMapping localToGlobalMapping;
  ierr = ISLocalToGlobalMappingCreate(mpiCommunicator(), 1, localSize(), 
                                      globalDofNos.data(), PETSC_COPY_VALUES, &localToGlobalMapping); CHKERRABORT(mpiCommunicator(),ierr);

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
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
globalSize()
{
  return globalSize_;
}
  
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
extractLocalNodes(std::vector<T> &vector)
{
  
}
  
template<int D, typename BasisFunctionType>
void MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
extractLocalDofs(std::vector<double> &vector)
{
  
}
  
}  // namespace