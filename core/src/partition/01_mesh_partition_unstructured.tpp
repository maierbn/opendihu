#include "partition/01_mesh_partition.h"

namespace Partition
{

// this is not yet implemented in parallel, it should however work serially
 
template<int D, typename BasisFunctionType>
MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
MeshPartition(global_no_t nElementsGlobal, global_no_t nNodesGlobal, std::shared_ptr<RankSubset> rankSubset):
  MeshPartitionBase(rankSubset), nElementsGlobal_(nElementsGlobal), nNodesGlobal_(nNodesGlobal) 
{
  this->initializeLocalDofsVector(nNodesLocalWithGhosts());
}

//! get the local to global mapping for the current partition
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
localToGlobalMapping()
{
  int nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>::nDofsPerNode();
  
  PetscErrorCode ierr;
  std::vector<PetscInt> globalDofNos(nNodesLocalWithGhosts()*nDofsPerNode);
  std::iota(globalDofNos.begin(), globalDofNos.end(), 0);
  ISLocalToGlobalMapping localToGlobalMapping;
  ierr = ISLocalToGlobalMappingCreate(mpiCommunicator(), 1, nNodesLocalWithGhosts()*nDofsPerNode, 
                                      globalDofNos.data(), PETSC_COPY_VALUES, &localToGlobalMapping); CHKERRABORT(mpiCommunicator(),ierr);

  return localToGlobalMapping;
}

//! number of entries in the current partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nElementsLocal()
{
  return nElementsGlobal_;
}

//! number of entries in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nElementsGlobal()
{
  return nElementsGlobal_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesLocalWithGhosts()
{
  return nNodesGlobal_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
nNodesGlobal()
{
  return nNodesGlobal_;
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
  
template<int D, typename BasisFunctionType>
void MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
output(std::ostream &stream)
{
  stream << "MeshPartition<Unstructured>, size global: " << globalSize_ << ", local: " << localSize_;
} 

}  // namespace