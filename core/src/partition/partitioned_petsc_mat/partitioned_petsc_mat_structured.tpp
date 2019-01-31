#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

#include "partition/mesh_partition/01_mesh_partition.h"
#include "petscis.h"

//! constructor, create square sparse matrix
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition, int nComponents,
                    int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, meshPartition, name),
  nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat<structured> (square sparse matrix) with " << nComponents_ << " components from meshPartition " << meshPartition;

  MatType matrixType = MATAIJ;  // sparse matrix type
  createMatrix(matrixType, diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor, create square sparse matrix
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition, int nComponents,
                    std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, meshPartition, name),
  nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat<structured> (square dense matrix) with " << nComponents_ << " components from meshPartition " << meshPartition;

  MatType matrixType = MATDENSE;  // dense matrix type
  createMatrix(matrixType, 0, 0);
}

//! constructor, create non-square matrix
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartitionRows,
                    std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns,
                    int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType>(meshPartitionRows, meshPartitionColumns, name),
  nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat<structured> (non-square sparse matrix) with " << nComponents_ << " components "
    << "from meshPartition rows: " << meshPartitionRows << ", columns: " << meshPartitionColumns;

  MatType matrixType = MATAIJ;  // sparse matrix type
  createMatrix(matrixType, diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor, create non-square matrix
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartitionRows,
                    std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns,
                    int nComponents, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType>(meshPartitionRows, meshPartitionColumns, name),
  nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat<structured> (non-square dense matrix) with " << nComponents_ << " components "
    << "from meshPartition rows: " << meshPartitionRows << ", columns: " << meshPartitionColumns;

  MatType matrixType = MATDENSE;  // dense matrix type
  createMatrix(matrixType, 0, 0);
}

//! constructor, use provided global matrix
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                    Mat &globalMatrix, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, meshPartition, name),
  nComponents_(1)
{
  VLOG(1) << "create PartitionedPetscMat<structured> from existing matrix, " << nComponents_ << " components from meshPartition " << meshPartition;

  assert(this->meshPartitionRows_);
  assert(this->meshPartitionColumns_);
  this->globalMatrix_ = globalMatrix;

  createLocalMatrix();
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
~PartitionedPetscMat()
{
  LOG(DEBUG) << "destroy PartitionedPetscMat";
  PetscErrorCode ierr;
  ierr = MatDestroy(&this->globalMatrix_); CHKERRV(ierr);
  ierr = MatDestroy(&this->localMatrix_); CHKERRV(ierr);
}

//! create a distributed Petsc matrix, according to the given partition
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
createMatrix(MatType matrixType, int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  dof_no_t nRowDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  dof_no_t nColumnDofsPerNode = ColumnsFunctionSpaceType::nDofsPerNode();

  dof_no_t nRowsLocal = this->meshPartitionRows_->nNodesLocalWithoutGhosts() * nRowDofsPerNode * nComponents_;
  dof_no_t nRowsGlobal = this->meshPartitionRows_->nNodesGlobal() * nRowDofsPerNode * nComponents_;

  dof_no_t nColumnsLocal = this->meshPartitionColumns_->nNodesLocalWithoutGhosts() * nColumnDofsPerNode * nComponents_;
  dof_no_t nColumnsGlobal = this->meshPartitionColumns_->nNodesGlobal() * nColumnDofsPerNode * nComponents_;

  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
  //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);

  assert(this->meshPartitionRows_);
  assert(this->meshPartitionColumns_);

  // parallel API
  //ierr = DMSetMatrixPreallocateOnly(this->dm_, PETSC_TRUE); CHKERRV(ierr);  // do not fill zero entries when DMCreateMatrix is called
  //ierr = DMCreateMatrix(dm_, &matrix_); CHKERRV(ierr);
  
  ierr = MatCreate(this->meshPartitionRows_->mpiCommunicator(), &this->globalMatrix_); CHKERRV(ierr);
  ierr = MatSetSizes(this->globalMatrix_, nRowsLocal, nColumnsLocal,
                      nRowsGlobal, nColumnsGlobal); CHKERRV(ierr);

  ierr = MatSetType(this->globalMatrix_, matrixType); CHKERRV(ierr);

  //ierr = MatSetFromOptions(this->globalMatrix_); CHKERRV(ierr);         // either use MatSetFromOptions or MatSetUp to allocate internal data structures
  ierr = PetscObjectSetName((PetscObject) this->globalMatrix_, this->name_.c_str()); CHKERRV(ierr);
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  // MATDENSE = "dense" - A matrix type to be used for dense matrices. This matrix type is identical to MATSEQDENSE when constructed with a single process communicator, and MATMPIDENSE otherwise.
  if (std::string(matrixType) == MATMPIDENSE || std::string(matrixType) == MATSEQDENSE || std::string(matrixType) == MATDENSE)
  {
    ierr = MatSetUp(this->globalMatrix_); CHKERRV(ierr);
  }
  else
  {
    // sparse matrix: preallocation of internal data structure
    // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
    // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
    // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
    // It is recommended that you call both of the above preallocation routines for simplicity.
    ierr = MatSeqAIJSetPreallocation(this->globalMatrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
    ierr = MatMPIAIJSetPreallocation(this->globalMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
    LOG(DEBUG) << "Mat SetPreallocation, diagonalNonZeros: " << diagonalNonZeros << ", offdiagonalNonZeros: " << offdiagonalNonZeros;
  }

  createLocalMatrix();
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
createLocalMatrix()
{
  VLOG(1) << "set local to global mapping for matrix: rows: " << (ISLocalToGlobalMapping)this->meshPartitionRows_->localToGlobalMappingDofs()
   << ", columns: " << (ISLocalToGlobalMapping)this->meshPartitionColumns_->localToGlobalMappingDofs();

  PetscErrorCode ierr;
  ierr = MatSetLocalToGlobalMapping(this->globalMatrix_, this->meshPartitionRows_->localToGlobalMappingDofs(), this->meshPartitionColumns_->localToGlobalMappingDofs()); CHKERRV(ierr);

  // output size of created matrix
  int nRows, nRowsLocal, nColumns, nColumnsLocal;
  ierr = MatGetSize(this->globalMatrix_, &nRows, &nColumns); CHKERRV(ierr);
  ierr = MatGetLocalSize(this->globalMatrix_, &nRowsLocal, &nColumnsLocal); CHKERRV(ierr);
  LOG(DEBUG) << "matrix \"" << this->name_ << "\" created, size global: " << nRows << "x" << nColumns << ", local: " << nRowsLocal << "x" << nColumnsLocal;

  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartitionRows_->dofNosLocalIS(), this->meshPartitionColumns_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" setValue " << (mode==INSERT_VALUES? "(insert)" : "(add)") << ", row " << row << ", col " << col << ", value " << value;
    VLOG(2) << stream.str();
  }

  // this wraps the standard PETSc MatSetValue on the local matrix
  PetscErrorCode ierr;
  ierr = MatSetValuesLocal(this->localMatrix_, 1, &row, 1, &col, &value, mode); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" setValues " << (addv==INSERT_VALUES? "(insert)" : "(add)") << ", rows [";
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
  ierr = MatSetValuesLocal(this->localMatrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
zeroRowsColumns(PetscInt numRows, const PetscInt rows[], PetscScalar diag)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" zeroRowsColumns rows [";
    for (int i = 0; i < numRows; i++)
      stream << rows[i] << " ";
    stream << "], diag " << diag;
    VLOG(2) << stream.str();
  }

  if (numRows == 0)
    return;
  
  PetscErrorCode ierr;
  // execute zeroRowsColumns on the global matrix, because it is not defined on the local matrix
  ierr = MatZeroRowsColumnsLocal(this->globalMatrix_, numRows, rows, diag, NULL, NULL); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartitionRows_->dofNosLocalIS(), this->meshPartitionColumns_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);

}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
zeroEntries()
{
  PetscErrorCode ierr;
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  ierr = MatZeroEntries(this->globalMatrix_); CHKERRV(ierr);
  
  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartitionRows_->dofNosLocalIS(), this->meshPartitionColumns_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
assembly(MatAssemblyType type)
{
  VLOG(2) << "\"" << this->name_ << "\" assembly " << (type == MAT_FINAL_ASSEMBLY? "(MAT_FINAL_ASSEMBLY)" : "(MAT_FLUSH_ASSEMBLY)");
  // type is either MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  
  // this wraps the standard PETSc assembleBegin/End
  PetscErrorCode ierr;
  
  // map the values of the localMatrix back into the global matrix
  ierr = MatRestoreLocalSubMatrix(this->globalMatrix_, this->meshPartitionRows_->dofNosLocalIS(), this->meshPartitionColumns_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, type); CHKERRV(ierr);
  
  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartitionRows_->dofNosLocalIS(), this->meshPartitionColumns_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;
  
  // map the values of the localMatrix back into the global matrix
  //ierr = MatRestoreLocalSubMatrix(this->globalMatrix_, this->meshPartitionRows_->dofNosLocalIS(),
  //                                this->meshPartitionColumns_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
  
  // assemble the global matrix
  //ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  //ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  if (VLOG_IS_ON(1))
  {
    for (int i = 0; i < m; i++)
    {
      assert(idxm[i] < this->meshPartitionRows_->nDofsGlobal());
    }
    for (int i = 0; i < n; i++)
    {
      assert(idxn[i] < this->meshPartitionColumns_->nDofsGlobal());
    }
  }

  // access the global matrix
  ierr = MatGetValues(this->globalMatrix_, m, idxm, n, idxn, v); CHKERRV(ierr);

  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" getValuesGlobalPetscIndexing, rows [";
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
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;

  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" getValues, rows [";
    for (int i = 0; i < m; i++)
      stream << idxm[i] << " ";
    stream << "], cols [";
    for (int i = 0; i < n; i++)
      stream << idxn[i] << " ";
    stream << "], ";
    VLOG(2) << stream.str();
  }
  // transfer the local indices to global indices
  std::vector<int> rowIndicesGlobal(m);
  std::vector<int> columnIndicesGlobal(n);
  ierr = ISLocalToGlobalMappingApply(this->meshPartitionRows_->localToGlobalMappingDofs(), m, idxm, rowIndicesGlobal.data()); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingApply(this->meshPartitionColumns_->localToGlobalMappingDofs(), n, idxn, columnIndicesGlobal.data()); CHKERRV(ierr);

  if (VLOG_IS_ON(2))
  {
    VLOG(2) << "rowIndicesGlobal: " << rowIndicesGlobal << ", columnIndicesGlobal: " << columnIndicesGlobal;
  }

  // access the global matrix
  ierr = MatGetValues(this->globalMatrix_, m, rowIndicesGlobal.data(), n, columnIndicesGlobal.data(), v); CHKERRV(ierr);

  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "values [";
    for (int i = 0; i < n*m; i++)
      stream << v[i] << " ";
    stream << "]";
    VLOG(2) << stream.str();
  }

}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
Mat &PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
valuesLocal()
{
  return this->localMatrix_;
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
Mat &PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
valuesGlobal()
{
  return this->globalMatrix_;
}

template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType,Mesh::isStructured<MeshType>>::
output(std::ostream &stream) const
{
#ifndef NDEBUG  
  // this method gets all values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo, nRanks;
  MPIUtility::handleReturnValue(MPI_Comm_rank(this->meshPartitionRows_->mpiCommunicator(), &ownRankNo), "MPI_Comm_rank");
  MPIUtility::handleReturnValue(MPI_Comm_size(this->meshPartitionRows_->mpiCommunicator(), &nRanks), "MPI_Comm_size");
  
  // get global size of matrix 
  int nRowsGlobal, nColumnsGlobal, nRowsLocal, nColumnsLocal;
  PetscErrorCode ierr;
  ierr = MatGetSize(this->globalMatrix_, &nRowsGlobal, &nColumnsGlobal); CHKERRV(ierr);
  ierr = MatGetLocalSize(this->globalMatrix_, &nRowsLocal, &nColumnsLocal); CHKERRV(ierr);
  
  // retrieve local values
  int nRowDofsLocal = this->meshPartitionRows_->nDofsLocalWithoutGhosts();
  int nColumnDofsLocal = this->meshPartitionColumns_->nDofsLocalWithoutGhosts();
  std::vector<double> localValues;
  localValues.resize(nRowDofsLocal*nColumnDofsLocal);
  
  getValues(nRowDofsLocal, this->meshPartitionRows_->dofNosLocal().data(), nColumnDofsLocal, this->meshPartitionColumns_->dofNosLocal().data(),
            localValues.data());
  
  // on every rank prepare a string with the local information
  std::string str;
  std::stringstream s;
  s << "Rank " << ownRankNo << " (" << *this->meshPartitionRows_->rankSubset() << "," << *this->meshPartitionColumns_->rankSubset() << "): "
    << PetscUtility::getStringMatrix(localValues, nRowsLocal, nColumnsLocal, nRowsGlobal, nColumnsGlobal);
  str = s.str();
  
  VLOG(1) << str;

  // exchange the lengths of the local information
  std::vector<int> localSizes(nRanks);
  localSizes[ownRankNo] = str.length();
  if (ownRankNo == 0)
  {
    MPIUtility::handleReturnValue(MPI_Gather(MPI_IN_PLACE, 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, this->meshPartitionRows_->mpiCommunicator()), "MPI_Gather (1)");
  }
  else
  {
    MPIUtility::handleReturnValue(MPI_Gather(localSizes.data(), 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, this->meshPartitionRows_->mpiCommunicator()), "MPI_Gather (1)");
  }
  
  // determine the maximum length
  int maxLocalSize;
  MPIUtility::handleReturnValue(MPI_Allreduce(localSizes.data() + ownRankNo, &maxLocalSize, 1, MPI_INT, MPI_MAX, this->meshPartitionRows_->mpiCommunicator()), "MPI_Allreduce");
  
  // send all string representations from the ranks to rank 0
  std::vector<char> recvBuffer(maxLocalSize*nRanks,0);
  std::vector<char> sendBuffer(maxLocalSize,0);
  memcpy(sendBuffer.data(), str.c_str(), str.length());
  
  // MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
  // Note that the recvcount argument at the root indicates the number of items it receives from each process, not the total number of items it receives.
  VLOG(1) << "MPI_Gather (2)";
  MPIUtility::handleReturnValue(MPI_Gather(sendBuffer.data(), maxLocalSize, MPI_CHAR, recvBuffer.data(), maxLocalSize, MPI_CHAR, 0, this->meshPartitionRows_->mpiCommunicator()), "MPI_Gather (2)");
  
  if (ownRankNo == 0)
  {
    // on rank 0 concatenate the strings from the ranks
    stream << "matrix \"" << this->name_ << "\" (" << nRowsGlobal << "x" << nColumnsGlobal << ")" << std::endl;
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      std::string s(recvBuffer.data() + maxLocalSize*rankNo, localSizes[rankNo]);
      stream << s << std::endl;
    }
  }

  PetscViewer viewer;
  static int counter = 0;
  std::stringstream matrixOutputFilename;
  matrixOutputFilename << "matrix_" << counter++ << ".txt";
  ierr = PetscViewerASCIIOpen(this->meshPartitionRows_->mpiCommunicator(), matrixOutputFilename.str().c_str(), &viewer); CHKERRV(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE); CHKERRV(ierr);
  ierr = MatView(this->globalMatrix_, viewer); CHKERRV(ierr);

  if (ownRankNo == 0)
  {
    stream << "(Dense matrix also written to \"" << matrixOutputFilename.str() << "\".)";
  }

  //ierr = MatView(this->globalMatrix_, PETSC_VIEWER_STDOUT_(this->meshPartitionRows_->mpiCommunicator())); CHKERRV(ierr);

#endif
}

template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<FunctionSpaceType> &matrix)
{
  matrix.output(stream);
  return stream;
}
