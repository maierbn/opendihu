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

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;
  
  // map the values of the localMatrix back into the global matrix
  ierr = MatRestoreLocalSubMatrix(this->globalMatrix_, this->meshPartition_->dofNosLocalIS(), this->meshPartition_->dofNosLocalIS(), &this->localMatrix_); CHKERRV(ierr);
  
  // assemble the global matrix
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  // access the global matrix
  ierr = MatGetValues(this->globalMatrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}

template<typename MeshType, typename BasisFunctionType>
Mat &PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
valuesLocal()
{
  return this->localMatrix_;
}

template<typename MeshType, typename BasisFunctionType>
Mat &PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
valuesGlobal()
{
  return this->globalMatrix_;
}

template<typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
output(std::ostream &stream) const
{
#ifndef NDEBUG  
  // this method gets all values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
  MPI_Comm_size(PETSC_COMM_WORLD, &nRanks);
  
  // get global size of matrix 
  int nRows, nColumns, nRowsLocal, nColumnsLocal;
  PetscErrorCode ierr;
  ierr = MatGetSize(this->globalMatrix_, &nRows, &nColumns); CHKERRV(ierr);
  ierr = MatGetLocalSize(this->globalMatrix_, &nRowsLocal, &nColumnsLocal); CHKERRV(ierr);
  
  // retrieve local values
  //int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
  //std::vector<double> localValues;
  //localValues.resize(MathUtility::sqr(nDofsLocal));
  
  //MatGetValuesLocal(this->globalMatrix_, nDofsLocal, this->meshPartition_->dofNosLocal(), nDofsLocal, this->meshPartition_->dofNosLocal(), localValues);
  
  // on every rank prepare a string with the local information
  std::string str;
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    if (rankNo == ownRankNo)
    {    
      std::stringstream s;
      s << "Rank " << rankNo << ": " << PetscUtility::getStringMatrix(this->globalMatrix_);
      str = s.str();
    }
  }
  
  // exchange the lengths of the local information
  std::vector<int> localSizes(nRanks);
  localSizes[ownRankNo] = str.length();
  MPI_Gather(localSizes.data() + ownRankNo, 1, MPI_INT, localSizes.data(), nRanks, MPI_INT, 0, this->meshPartition_->mpiCommunicator());
  
  // determine the maximum length
  int maxLocalSize;
  MPI_Allreduce(localSizes.data() + ownRankNo, &maxLocalSize, 1, MPI_INT, MPI_MAX, this->meshPartition_->mpiCommunicator());
  
  // send all string representations from the ranks to rank 0
  std::vector<char> recvBuffer(maxLocalSize*nRanks,0);
  std::vector<char> sendBuffer(maxLocalSize,0);
  memcpy(sendBuffer.data(), str.c_str(), str.length());
  
  MPI_Gather(sendBuffer.data(), maxLocalSize, MPI_CHAR, recvBuffer.data(), maxLocalSize*nRanks, MPI_CHAR, 0, this->meshPartition_->mpiCommunicator());
  
  if (ownRankNo == 0)
  {
    // on rank 0 concatenate the strings from the ranks
    stream << "matrix \"" << this->name_ << "\" (" << nRows << "x" << nColumns << ")" << std::endl;
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      std::string s(recvBuffer.data() + maxLocalSize*rankNo, localSizes[rankNo]);
      stream << s << std::endl;
    }
  }
#endif
}

template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<BasisOnMeshType> &matrix)
{
  matrix.output(stream);
  return stream;
}
