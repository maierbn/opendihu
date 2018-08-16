#include "partition/partitioned_petsc_mat.h"

//! constructor
template<int D, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>> meshPartition, 
                    int nComponents, int diagonalNonZeros, int offdiagonalNonZeros) :
  meshPartition_(meshPartition), nComponents_(nComponents)
{
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType> BasisOnMeshType;
  
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor
template<int D, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(Mat &matrix) :
  meshPartition_(nullptr), nComponents_(-1)
{ 
  //! constructor to simply wrap an existing Mat, as needed in nonlinear solver callback functions for jacobians
  this->matrix_ = matrix;
}

//! create a distributed Petsc matrix, according to the given partition
template<int D, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  assert(meshPartition_);
  
  dof_no_t nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>::nDofsPerNode();
  dof_no_t nRowsLocal = this->meshPartition_->nLocalNodes() * nDofsPerNode;
  dof_no_t nRowsGlobal = this->meshPartition_->nGlobalNodes() * nDofsPerNode;
  
  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
  //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);
  
  ierr = MatCreate(meshPartition_->mpiCommunicator(), &matrix_); CHKERRV(ierr);
  ierr = MatSetSizes(matrix_, nRowsLocal, nRowsLocal, nRowsGlobal, nRowsGlobal); CHKERRV(ierr);
  ierr = MatSetFromOptions(matrix_); CHKERRV(ierr);                        
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

  // sparse matrix
  ierr = MatMPIAIJSetPreallocation(matrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(matrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
 
  ierr = MatSetLocalToGlobalMapping(matrix_, meshPartition_->localToGlobalMapping(), meshPartition_->localToGlobalMapping()); CHKERRV(ierr);
}