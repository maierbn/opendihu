#include "partition/partitioned_petsc_mat.h"

#include "partition/01_mesh_partition.h"

//! constructor
template<typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>> meshPartition, int nComponents,
                    int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatBase<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>(meshPartition, name), nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat<structured> with " << nComponents_ << " components from meshPartition " << meshPartition;
  
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! create a distributed Petsc matrix, according to the given partition
template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  dof_no_t nDofsPerNode = BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>::nDofsPerNode();
  dof_no_t nRowsLocal = this->meshPartition_->nNodesLocalWithoutGhosts() * nDofsPerNode * nComponents_;
  dof_no_t nRowsGlobal = this->meshPartition_->nNodesGlobal() * nDofsPerNode * nComponents_;
  
  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
  //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);

  assert(this->meshPartition_);

  // parallel API
  //ierr = DMSetMatrixPreallocateOnly(this->dm_, PETSC_TRUE); CHKERRV(ierr);  // do not fill zero entries when DMCreateMatrix is called
  //ierr = DMCreateMatrix(dm_, &matrix_); CHKERRV(ierr);
  
  ierr = MatCreate(this->meshPartition_->mpiCommunicator(), &this->globalMatrix_); CHKERRV(ierr);
  ierr = MatSetSizes(this->globalMatrix_, nRowsLocal, nRowsLocal, 
                      nRowsGlobal, nRowsGlobal); CHKERRV(ierr);
  ierr = MatSetFromOptions(this->globalMatrix_); CHKERRV(ierr);         // either use MatSetFromOptions or MatSetUp to allocate internal data structures                  
  ierr = PetscObjectSetName((PetscObject) this->globalMatrix_, this->name_.c_str()); CHKERRV(ierr);
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

  // sparse matrix: preallocation of internal data structure
  ierr = MatMPIAIJSetPreallocation(this->globalMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(this->globalMatrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
  
  VLOG(1) << "set local to global mapping for matrix: " << (ISLocalToGlobalMapping)this->meshPartition_->localToGlobalMappingDofs();
  
  ierr = MatSetLocalToGlobalMapping(this->globalMatrix_, this->meshPartition_->localToGlobalMappingDofs(), this->meshPartition_->localToGlobalMappingDofs()); CHKERRV(ierr);
  
  // output size of create matrix
  int nRows, nColumns, nColumnsLocal;
  ierr = MatGetSize(this->globalMatrix_, &nRows, &nColumns); CHKERRV(ierr);
  ierr = MatGetLocalSize(this->globalMatrix_, &nRowsLocal, &nColumnsLocal); CHKERRV(ierr);
  LOG(DEBUG) << "matrix created, size global: " << nRows << "x" << nColumns << ", local: " << nRowsLocal << "x" << nColumnsLocal;

  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<BasisOnMeshType> &matrix)
{
  matrix.output(stream);
  return stream;
}
