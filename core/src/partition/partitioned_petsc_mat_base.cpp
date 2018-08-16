#include "partition/partitioned_petsc_mat_base.h"

void PartitionedPetscMatBase::
setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
{
  // this wraps the standard PETSc MatSetValue on the local matrix
  PetscErrorCode ierr;
  
  ierr = MatSetValuesLocal(matrix_, 1, &row, 1, &col, &value, mode); CHKERRV(ierr);
}

void PartitionedPetscMatBase::
setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  // this wraps the standard PETSc MatSetValues on the local matrix
  PetscErrorCode ierr;
  
  ierr = MatSetValuesLocal(matrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

void PartitionedPetscMatBase::
zeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag)
{
  PetscErrorCode ierr;
  ierr = MatZeroRowsColumnsLocal(matrix_, numRows, rows, diag, NULL, NULL); CHKERRV(ierr);
}

void PartitionedPetscMatBase::
zeroEntries()
{
  PetscErrorCode ierr;
  ierr = MatZeroEntries(matrix_); CHKERRV(ierr);  
}

void PartitionedPetscMatBase::
assembly(MatAssemblyType type)
{
  // this wraps the standard PETSc assembleBegin/End
  PetscErrorCode ierr;
  
  ierr = MatAssemblyBegin(matrix_, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrix_, type); CHKERRV(ierr);
}

void PartitionedPetscMatBase::
getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;
  
  ierr = MatGetValues(matrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

//! get a reference to the PETSc matrix
Mat &PartitionedPetscMatBase::
values()
{
  return this->matrix_;
}