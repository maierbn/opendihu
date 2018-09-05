#include "partition/partitioned_petsc_mat.h"

#include "partition/01_mesh_partition.h"

//! constructor
template<typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition, int nComponents,
                    int diagonalNonZeros, int offdiagonalNonZeros, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, name), nComponents_(nComponents)
{
  VLOG(1) << "create PartitionedPetscMat<structured> with " << nComponents_ << " components from meshPartition " << meshPartition;

  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor, use provided global matrix
template<typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                    Mat &globalMatrix, std::string name) :
  PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, name), nComponents_(1)
{
  VLOG(1) << "create PartitionedPetscMat<structured> from existing matrix, " << nComponents_ << " components from meshPartition " << meshPartition;

  assert(this->meshPartition_);
  this->globalMatrix_ = globalMatrix;

  createLocalMatrix();
}

//! create a distributed Petsc matrix, according to the given partition
template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  dof_no_t nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
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

  createLocalMatrix();
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
createLocalMatrix()
{
  VLOG(1) << "set local to global mapping for matrix: " << (ISLocalToGlobalMapping)this->meshPartition_->localToGlobalMappingDofs();

  PetscErrorCode ierr;
  ierr = MatSetLocalToGlobalMapping(this->globalMatrix_, this->meshPartition_->localToGlobalMappingDofs(), this->meshPartition_->localToGlobalMappingDofs()); CHKERRV(ierr);

  // output size of create matrix
  int nRows, nRowsLocal, nColumns, nColumnsLocal;
  ierr = MatGetSize(this->globalMatrix_, &nRows, &nColumns); CHKERRV(ierr);
  ierr = MatGetLocalSize(this->globalMatrix_, &nRowsLocal, &nColumnsLocal); CHKERRV(ierr);
  LOG(DEBUG) << "matrix created, size global: " << nRows << "x" << nColumns << ", local: " << nRowsLocal << "x" << nColumnsLocal;

  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
    VLOG(2) << stream.str();
  }
  
  // this wraps the standard PETSc MatSetValues on the local matrix
  PetscErrorCode ierr;
  ierr = MatSetValuesLocal(this->localMatrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
  
  PetscErrorCode ierr;
  // execute zeroRowsColumns on the global matrix, because it is not defined on the local matrix
  ierr = MatZeroRowsColumnsLocal(this->globalMatrix_, numRows, rows, diag, NULL, NULL); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);

}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
zeroEntries()
{
  PetscErrorCode ierr;
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  
  ierr = MatZeroEntries(this->globalMatrix_); CHKERRV(ierr);
  
  // get the local submatrix from the global matrix
  ierr = MatGetLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
assembly(MatAssemblyType type)
{
  VLOG(2) << "\"" << this->name_ << "\" assembly " << (type == MAT_FINAL_ASSEMBLY? "(MAT_FINAL_ASSEMBLY)" : "(MAT_FLUSH_ASSEMBLY)");
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

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;
  
  // map the values of the localMatrix back into the global matrix
  ierr = MatRestoreLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(),
                                  this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  // access the global matrix
  ierr = MatGetValues(this->globalMatrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;

  // transfer the local indices to global indices
  std::vector<int> rowIndicesGlobal(m);
  std::vector<int> columnIndicesGlobal(n);
  ierr = ISLocalToGlobalMappingApply(this->meshPartition_->localToGlobalMappingDofs(), m, idxm, rowIndicesGlobal.data()); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingApply(this->meshPartition_->localToGlobalMappingDofs(), n, idxn, columnIndicesGlobal.data()); CHKERRV(ierr);

  // access the global matrix
  ierr = MatGetValues(this->globalMatrix_, m, rowIndicesGlobal.data(), n, columnIndicesGlobal.data(), v); CHKERRV(ierr);

}

template<typename MeshType, typename BasisFunctionType>
Mat &PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
valuesLocal()
{
  return this->localMatrix_;
}

template<typename MeshType, typename BasisFunctionType>
Mat &PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
valuesGlobal()
{
  return this->globalMatrix_;
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
output(std::ostream &stream) const
{
#ifndef NDEBUG  
  // this method gets all values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo, nRanks;
  MPI_Comm_rank(this->meshPartition_->mpiCommunicator(), &ownRankNo);
  MPI_Comm_size(this->meshPartition_->mpiCommunicator(), &nRanks);
  
  // get global size of matrix 
  int nRowsGlobal, nColumnsGlobal, nRowsLocal, nColumnsLocal;
  PetscErrorCode ierr;
  ierr = MatGetSize(this->globalMatrix_, &nRowsGlobal, &nColumnsGlobal); CHKERRV(ierr);
  ierr = MatGetLocalSize(this->globalMatrix_, &nRowsLocal, &nColumnsLocal); CHKERRV(ierr);
  
  // retrieve local values
  int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
  std::vector<double> localValues;
  localValues.resize(MathUtility::sqr(nDofsLocal));
  
  getValues(nDofsLocal, this->meshPartition_->dofNosLocal().data(), nDofsLocal, this->meshPartition_->dofNosLocal().data(),
            localValues.data());

  
  // on every rank prepare a string with the local information
  std::string str;
  std::stringstream s;
  s << "Rank " << ownRankNo << ": " << PetscUtility::getStringMatrix(localValues, nRowsLocal, nColumnsLocal, nRowsGlobal, nColumnsGlobal);
  str = s.str();
  
  //VLOG(1) << str;


  // exchange the lengths of the local information
  std::vector<int> localSizes(nRanks);
  localSizes[ownRankNo] = str.length();
  MPI_Gather(localSizes.data() + ownRankNo, 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, this->meshPartition_->mpiCommunicator());
  
  // determine the maximum length
  int maxLocalSize;
  MPI_Allreduce(localSizes.data() + ownRankNo, &maxLocalSize, 1, MPI_INT, MPI_MAX, this->meshPartition_->mpiCommunicator());
  
  // send all string representations from the ranks to rank 0
  std::vector<char> recvBuffer(maxLocalSize*nRanks,0);
  std::vector<char> sendBuffer(maxLocalSize,0);
  memcpy(sendBuffer.data(), str.c_str(), str.length());
  
  // MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
  // Note that the recvcount argument at the root indicates the number of items it receives from each process, not the total number of items it receives.
  MPI_Gather(sendBuffer.data(), maxLocalSize, MPI_CHAR, recvBuffer.data(), maxLocalSize, MPI_CHAR, 0, this->meshPartition_->mpiCommunicator());
  
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
  ierr = PetscViewerASCIIOpen(this->meshPartition_->mpiCommunicator(), matrixOutputFilename.str().c_str(), &viewer); CHKERRV(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE); CHKERRV(ierr);
  ierr = MatView(this->globalMatrix_, viewer); CHKERRV(ierr);

  if (ownRankNo == 0)
  {
    stream << "(Dense matrix also written to \"" << matrixOutputFilename.str() << "\".)";
  }

  //ierr = MatView(this->globalMatrix_, PETSC_VIEWER_STDOUT_(this->meshPartition_->mpiCommunicator())); CHKERRV(ierr);

#endif
}

template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<FunctionSpaceType> &matrix)
{
  matrix.output(stream);
  return stream;
}
