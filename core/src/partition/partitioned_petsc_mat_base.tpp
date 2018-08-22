#include "partition/partitioned_petsc_mat_base.h"


template<typename BasisOnMeshType>
PartitionedPetscMatBase<BasisOnMeshType>::
PartitionedPetscMatBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition, std::string name) :
  meshPartition_(meshPartition), name_(name)
{
}
  
template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "matrix \"" << name_ << "\" setValue " << (mode==INSERT_VALUES? "(insert)" : "(add)") << ", row " << row << ", col " << col << ", value " << value;
    VLOG(2) << stream.str();
  }
  
  // this wraps the standard PETSc MatSetValue on the local matrix
  PetscErrorCode ierr;
  ierr = MatSetValuesLocal(this->localMatrix_, 1, &row, 1, &col, &value, mode); CHKERRV(ierr);
}

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "matrix \"" << name_ << "\" setValues " << (addv==INSERT_VALUES? "(insert)" : "(add)") << ", rows [";
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

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
zeroRowsColumns(PetscInt numRows, const PetscInt rows[], PetscScalar diag)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "matrix \"" << name_ << "\" zeroRowsColumns rows [";
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
  VLOG(2) << "matrix \"" << name_ << "\" assembly " << (type == MAT_FINAL_ASSEMBLY? "(MAT_FINAL_ASSEMBLY)" : "(MAT_FLUSH_ASSEMBLY)");
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
  ierr = MatAssemblyBegin(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->globalMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
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

template<typename BasisOnMeshType>
void PartitionedPetscMatBase<BasisOnMeshType>::
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