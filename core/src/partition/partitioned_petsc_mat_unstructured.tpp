#include "partition/partitioned_petsc_mat.h"

//! constructor
template<int D, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>> meshPartition, 
                    int nComponents, int diagonalNonZeros, int offdiagonalNonZeros) :
  PartitionedPetscMatBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>(meshPartition), nComponents_(nComponents)
{
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType> BasisOnMeshType;
  
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! create a distributed Petsc matrix, according to the given partition
template<int D, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  assert(this->meshPartition_);
  
  dof_no_t nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>::nDofsPerNode();
  dof_no_t nRowsLocal = this->meshPartition_->nNodesLocalWithoutGhosts() * nDofsPerNode;
  dof_no_t nRowsGlobal = this->meshPartition_->nNodesGlobal() * nDofsPerNode;
  
  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
  //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);
  
  ierr = MatCreate(this->meshPartition_->mpiCommunicator(), &this->globalMatrix_); CHKERRV(ierr);
  ierr = MatSetSizes(this->globalMatrix_, nRowsLocal, nRowsLocal, nRowsGlobal, nRowsGlobal); CHKERRV(ierr);
  ierr = MatSetFromOptions(this->globalMatrix_); CHKERRV(ierr);                        
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

  // sparse matrix
  ierr = MatMPIAIJSetPreallocation(this->globalMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(this->globalMatrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
 
  ierr = MatSetLocalToGlobalMapping(this->globalMatrix_, this->meshPartition_->localToGlobalMapping(), this->meshPartition_->localToGlobalMapping()); CHKERRV(ierr);
}