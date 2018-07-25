#include "partition/partitioned_petsc_vec.h"

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(std::shared_ptr<MeshPartition> meshPartition, int nComponents,
                    int diagonalNonZeros, int offdiagonalNonZeros) :
  meshPartition_(meshPartition), nComponents_(nComponents)
{
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType> BasisOnMeshType;
  
  PetscErrorCode ierr;
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(Mat &matrix) :
  meshPartition_(nullptr), nComponents_(-1), matrix_(matrix)
{ 
  //! constructor to simply wrap an existing Mat, as needed in nonlinear solver callback functions for jacobians
}

//! create a distributed Petsc matrix, according to the given partition
template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.localSize(), partition.localSize(), n, n,
  //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);

  assert(meshPartition_);
  
  ierr = MatCreate(meshPartition_->mpiCommunicator(), &matrix_); CHKERRV(ierr);
  ierr = MatSetSizes(matrix_, meshPartition_->localSize(), meshPartition_->localSize(), meshPartition->globalSize(), meshPartition->globalSize()); CHKERRV(ierr);
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