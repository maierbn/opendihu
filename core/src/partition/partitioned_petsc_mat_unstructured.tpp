#include "partition/partitioned_petsc_mat.h"

//! constructor, create square matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>> meshPartition,
                    int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartition, meshPartition, name),
  nComponents_(nComponents)
{
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor, create non-square matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionRows,
                    std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionColumns,
                    int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartitionRows, meshPartitionColumns, name),
  nComponents_(nComponents)
{
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}


//! constructor, use provided global matrix
template<int D, typename BasisFunctionType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                      Mat &globalMatrix, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(meshPartition, meshPartition, name),
  nComponents_(1)
{
  VLOG(1) << "create PartitionedPetscMat<unstructured> from existing matrix, " << nComponents_ << " components from meshPartition " << meshPartition;

  assert(this->meshPartitionRows_);
  assert(this->meshPartitionColumns_);
  this->matrix_ = globalMatrix;
}

//! create a distributed Petsc matrix, according to the given partition
template<int D, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
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
  ierr = MatSetFromOptions(this->matrix_); CHKERRV(ierr);                        
  ierr = PetscObjectSetName((PetscObject) this->matrix_, this->name_.c_str()); CHKERRV(ierr);
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

  // sparse matrix
  ierr = MatMPIAIJSetPreallocation(this->matrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(this->matrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
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
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
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
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
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
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
zeroEntries()
{
  PetscErrorCode ierr;
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->matrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  ierr = MatZeroEntries(this->matrix_); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
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
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
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
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;

  // access the global matrix
  ierr = MatGetValues(this->matrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

template<int D, typename BasisFunctionType>
Mat &PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
valuesLocal()
{
  return this->matrix_;
}

template<int D, typename BasisFunctionType>
Mat &PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
valuesGlobal()
{
  return this->matrix_;
}

template<int D, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>::
output(std::ostream &stream) const
{
  // this method gets all values and outputs them to stream
  PetscMPIInt ownRankNo, nRanks;
  MPI_Comm_rank(this->meshPartitionRows_->mpiCommunicator(), &ownRankNo);
  MPI_Comm_size(this->meshPartitionRows_->mpiCommunicator(), &nRanks);
  
  // get global size of matrix 
  int nRows, nColumns;
  PetscErrorCode ierr;
  ierr = MatGetSize(this->matrix_, &nRows, &nColumns); CHKERRV(ierr);
  
  ierr = MatAssemblyBegin(this->matrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->matrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  stream << PetscUtility::getStringMatrix(this->matrix_) << std::endl;
}
