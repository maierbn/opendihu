#include "partition/partitioned_petsc_mat.h"

#include "partition/01_mesh_partition.h"

//! constructor
template<typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>> meshPartition, int nComponents,
                    int diagonalNonZeros, int offdiagonalNonZeros) :
  meshPartition_(meshPartition), nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat from meshPartition " << meshPartition << " with " << nComponents_ << " components";
  
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor
template<typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(Mat &matrix) :
  meshPartition_(nullptr), nComponents_(-1)
{ 
  VLOG(1) << "create PartitionedPetscMat as wrapper to a Petsc Mat";
  
  //! constructor to simply wrap an existing Mat, as needed in nonlinear solver callback functions for jacobians
  this->matrix_ = matrix;
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
  
  const bool serial = true;   /// use the serial PETSc API
  
  // sparse matrix, serial
  if (serial)
  {
    //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
    //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);

    assert(meshPartition_);
   
    ierr = MatCreate(meshPartition_->mpiCommunicator(), &matrix_); CHKERRV(ierr);
    ierr = MatSetSizes(matrix_, nRowsLocal, nRowsLocal, 
                       nRowsGlobal, nRowsGlobal); CHKERRV(ierr);
    ierr = MatSetFromOptions(matrix_); CHKERRV(ierr);                        
    
    // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
    //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

    // dense matrix
    //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

    // sparse matrix
    ierr = MatMPIAIJSetPreallocation(matrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
    ierr = MatSeqAIJSetPreallocation(matrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
    
    VLOG(1) << "set local to global mapping for matrix: " << (ISLocalToGlobalMapping)meshPartition_->localToGlobalMappingDofs();
    
    ierr = MatSetLocalToGlobalMapping(matrix_, meshPartition_->localToGlobalMappingDofs(), meshPartition_->localToGlobalMappingDofs()); CHKERRV(ierr);
  }
  else 
  {
    // parallel API
    //ierr = DMSetMatrixPreallocateOnly(this->dm_, PETSC_TRUE); CHKERRV(ierr);  // do not fill zero entries when DMCreateMatrix is called
    //ierr = DMCreateMatrix(dm_, &matrix_); CHKERRV(ierr);
  }
}