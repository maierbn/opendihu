#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component.h"

//! constructor, create square sparse matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>> meshPartition,
                                int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartition, meshPartition, name)
{
  MatType matrixType = MATAIJ;  // sparse matrix type
  createMatrix(matrixType, diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor, create square dense matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>> meshPartition,
                                std::string name) :
  PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartition, meshPartition, name)
{
  MatType matrixType = MATDENSE;  // dense matrix type
  createMatrix(matrixType, 0, 0);
}

//! constructor, create non-square sparse matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionRows,
                                std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionColumns,
                                int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartitionRows, meshPartitionColumns, name)
{
  MatType matrixType = MATAIJ;  // sparse matrix type
  createMatrix(matrixType, diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor, create non-square dense matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionRows,
                                std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionColumns,
                                std::string name) :
  PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartitionRows, meshPartitionColumns, name)
{
  MatType matrixType = MATDENSE;  // dense matrix type
  createMatrix(matrixType, 0, 0);
}

//! constructor, use provided global matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                                Mat &globalMatrix, std::string name) :
  PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartition, meshPartition, name)
{
  VLOG(1) << "create PartitionedPetscMatOneComponent<unstructured> from existing matrix, from meshPartition " << meshPartition;

  assert(this->meshPartitionRows_);
  assert(this->meshPartitionColumns_);
  this->matrix_ = globalMatrix;
}

//! create a distributed Petsc matrix, according to the given partition
template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
createMatrix(MatType matrixType, int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  assert(this->meshPartitionRows_);
  assert(this->meshPartitionColumns_);
  
  dof_no_t nRowDofs = this->meshPartitionRows_->nDofs();
  dof_no_t nColumnDofs = this->meshPartitionColumns_->nDofs();
  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
  //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);
  
  ierr = MatCreate(this->meshPartitionRows_->mpiCommunicator(), &this->matrix_); CHKERRV(ierr);
  ierr = MatSetSizes(this->matrix_, nRowDofs, nColumnDofs, nRowDofs, nColumnDofs); CHKERRV(ierr);

  ierr = MatSetType(this->matrix_, matrixType); CHKERRV(ierr);

  //ierr = MatSetFromOptions(this->matrix_); CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) this->matrix_, this->name_.c_str()); CHKERRV(ierr);
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  // MATDENSE = "dense" - A matrix type to be used for dense matrices. This matrix type is identical to MATSEQDENSE when constructed with a single process communicator, and MATMPIDENSE otherwise.
  if (std::string(matrixType) == MATMPIDENSE || std::string(matrixType) == MATSEQDENSE || std::string(matrixType) == MATDENSE)
  {
    ierr = MatSetUp(this->matrix_); CHKERRV(ierr);
  }
  else
  {
    // sparse matrix: preallocation of internal data structure
    // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
    // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
    // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
    // It is recommended that you call both of the above preallocation routines for simplicity.
    ierr = MatMPIAIJSetPreallocation(this->matrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
    ierr = MatSeqAIJSetPreallocation(this->matrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
  }
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << " \"" << this->name_ << "\" setValue " << (mode==INSERT_VALUES? "(insert)" : "(add)") << ", row " << row << ", col " << col << ", value " << value;
    VLOG(2) << stream.str();
  }
  
  // this wraps the standard PETSc MatSetValue on the local matrix
  PetscErrorCode ierr;
  ierr = MatSetValues(this->matrix_, 1, &row, 1, &col, &value, mode); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << " \"" << this->name_ << "\" setValues " << (addv==INSERT_VALUES? "(insert)" : "(add)") << ", rows [";
    for (int i = 0; i < m; i++)
      stream << idxm[i] << " ";
    stream << "], cols [";
    for (int i = 0; i < n; i++)
      stream << idxn[i] << " ";
    stream << "], values [";
    for (int i = 0; i < n*m; i++)
      stream << v[i] << " ";
    stream << "]";
    VLOG(2) << stream.str();
  }
  
  // this wraps the standard PETSc MatSetValues on the local matrix
  PetscErrorCode ierr;
  ierr = MatSetValues(this->matrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
zeroRowsColumns(PetscInt numRows, const PetscInt rows[], PetscScalar diag)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << " \"" << this->name_ << "\" zeroRowsColumns rows [";
    for (int i = 0; i < numRows; i++)
      stream << rows[i] << " ";
    stream << "], diag " << diag;
    VLOG(2) << stream.str();
  }
  
  PetscErrorCode ierr;
  ierr = MatZeroRowsColumns(this->matrix_, numRows, rows, diag, NULL, NULL); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->matrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
zeroRows(PetscInt numRows, const PetscInt rows[])
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" zeroRows rows [";
    for (int i = 0; i < numRows; i++)
      stream << rows[i] << " ";
    stream << "]";
    VLOG(2) << stream.str();
  }

  PetscErrorCode ierr;

  if (numRows != 0)
  {
    // this call is collective, but apparently fails if numRows == 0
    ierr = MatZeroRowsLocal(this->matrix_, numRows, rows, 0.0, NULL, NULL); CHKERRV(ierr);
  }
  else
  {
    std::vector<PetscInt> rows = {0};
    ierr = MatZeroRowsLocal(this->matrix_, 0, rows.data(), 0.0, NULL, NULL); CHKERRV(ierr);
  }
}
template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
zeroEntries()
{
  PetscErrorCode ierr;
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->matrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  ierr = MatZeroEntries(this->matrix_); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
assembly(MatAssemblyType type)
{
  VLOG(2) << " \"" << this->name_ << "\" assembly " << (type == MAT_FINAL_ASSEMBLY? "(MAT_FINAL_ASSEMBLY)" : "(MAT_FLUSH_ASSEMBLY)");
  // type is either MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  
  // this wraps the standard PETSc assembleBegin/End
  PetscErrorCode ierr;
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->matrix_, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, type); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;

  // assemble the global matrix
  ierr = MatAssemblyBegin(this->matrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  // access the global matrix
  ierr = MatGetValues(this->matrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
setValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  // this wraps the standard PETSc MatSetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;

  // access the global matrix
  ierr = MatSetValues(this->matrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;

  // access the global matrix
  ierr = MatGetValues(this->matrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
Mat &PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
valuesLocal()
{
  return this->matrix_;
}

template<int D, typename BasisFunctionType>
Mat &PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
valuesGlobal()
{
  return this->matrix_;
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
dumpMatrix(std::string filename, std::string format)
{
  // std::string filename, std::string format, Mat &matrix, MPI_Comm mpiCommunicator
  PetscUtility::dumpMatrix(filename, format, this->matrix_, this->meshPartitionRows_->mpiCommunicator());
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>::
output(std::ostream &stream) const
{
  // this method gets all values and outputs them to stream
  // get global size of matrix 
  PetscInt nRows, nColumns;
  PetscErrorCode ierr;
  ierr = MatGetSize(this->matrix_, &nRows, &nColumns); CHKERRV(ierr);
  
  ierr = MatAssemblyBegin(this->matrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  stream << PetscUtility::getStringMatrix(this->matrix_) << std::endl;
}
