#include "partition/partitioned_petsc_mat_base.h"


template<typename BasisOnMeshType>
PartitionedPetscMatBase<BasisOnMeshType>::
PartitionedPetscMatBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition) :
  meshPartition_(meshPartition)
{
}
  
template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
{
  // this wraps the standard PETSc MatSetValue on the local matrix
  PetscErrorCode ierr;
  
  ierr = MatSetValuesLocal(this->localMatrix_, 1, &row, 1, &col, &value, mode); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  // this wraps the standard PETSc MatSetValues on the local matrix
  PetscErrorCode ierr;
  
  ierr = MatSetValuesLocal(this->localMatrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
zeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag)
{
  PetscErrorCode ierr;
  ierr = MatZeroRowsColumnsLocal(this->localMatrix_, numRows, rows, diag, NULL, NULL); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
zeroEntries()
{
  PetscErrorCode ierr;
  ierr = MatZeroEntries(this->localMatrix_); CHKERRV(ierr);  
}

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
assembly(MatAssemblyType type)
{
  // type is either MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  
  // this wraps the standard PETSc assembleBegin/End
  PetscErrorCode ierr;
  
  // map the values of the localMatrix back into the global matrix
  ierr = MatRestoreLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, type); CHKERRV(ierr);
  
  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;
  
  // map the values of the localMatrix back into the global matrix
  ierr = MatRestoreLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  // access the global matrix
  ierr = MatGetValues(this->globalMatrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
//! get a reference to the PETSc matrix
Mat &PartitionedPetscMatBase<BasisOnMeshType>::
valuesLocal()
{
  return this->localMatrix_;
}

template<typename BasisOnMeshType>
//! get a reference to the PETSc matrix
Mat &PartitionedPetscMatBase<BasisOnMeshType>::
valuesGlobal()
{
  return this->globalMatrix_;
}