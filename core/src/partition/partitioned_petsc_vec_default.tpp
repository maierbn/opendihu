#include "partition/partitioned_petsc_vec.h"

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
PartitionedPetscVec<FunctionSpaceType,nComponents,DummyForTraits>::
PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, std::string name) :
  PartitionedPetscVecBase<FunctionSpaceType>(meshPartition, name)
{
  createVector();
}

//! constructor from existing vector, values are not copied
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
PartitionedPetscVec<FunctionSpaceType,nComponents,DummyForTraits>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpaceType,nComponents> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpaceType>(rhs.meshPartition(), name)
 
{
  createVector();
  
  // copy existing values
  PetscErrorCode ierr;
  for (int i = 0; i < nComponents; i++)
  {
    ierr = VecCopy(rhs.valuesGlobal(i), values_[i]); CHKERRV(ierr);
  }
}

//! constructor from existing vector, values are not copied
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
template<int nComponents2>
PartitionedPetscVec<FunctionSpaceType,nComponents,DummyForTraits>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpaceType,nComponents2> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpaceType>(rhs.meshPartition(), name)
 
{
  createVector();
  
  // copy existing values
  PetscErrorCode ierr;
  for (int i = 0; i < std::min(nComponents,nComponents2); i++)
  {
    ierr = VecCopy(rhs.valuesGlobal(i), values_[i]); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
createVector()
{
  assert(this->meshPartition_);
  
  dof_no_t nEntriesLocal = this->meshPartition_->nDofs();
  dof_no_t nEntriesGlobal = nEntriesLocal;
  
  PetscErrorCode ierr;

  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // initialize PETSc vector object
    ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &values_[componentNo]); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject) values_[componentNo], this->name_.c_str()); CHKERRV(ierr);

    // initialize size of vector
    ierr = VecSetSizes(values_[componentNo], nEntriesLocal, nEntriesGlobal); CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(values_[componentNo]); CHKERRV(ierr);
  }
}

//! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
startGhostManipulation()
{
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
finishGhostManipulation()
{
  assert(values_.size() == nComponents);
  
  PetscErrorCode ierr;
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecAssemblyBegin(values_[componentNo]); CHKERRV(ierr); 
  }
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecAssemblyEnd(values_[componentNo]); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
zeroGhostBuffer()
{
}

//! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  assert(componentNo >= 0);
  assert(componentNo < values_.size());
  assert(values_.size() == nComponents);
  
  // debugging output
  if (VLOG_IS_ON(3))
  {
    std::stringstream str;
    str << "\"" << this->name_ << "\" setValues(componentNo=" << componentNo << ", indices=";
    for (int i = 0; i < ni; i++)
    {
      str << ix[i] << " ";
    }
    str << "): ";
    for (int i = 0; i < ni; i++)
    {
      str << y[i] << " ";
    }
    VLOG(3) << str.str();
  }
  assert(!valuesContiguousInUse_);
  
  // this wraps the standard PETSc VecSetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(values_[componentNo], ni, ix, y, iora); CHKERRV(ierr);
}

//! wrapper to the PETSc VecSetValue, acting only on the local data
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  assert(componentNo >= 0);
  assert(componentNo < values_.size());
  assert(values_.size() == nComponents);
  assert(!valuesContiguousInUse_);

  // debugging output
  VLOG(3) << "\"" << this->name_ << "\" setValue(componentNo=" << componentNo << ", row=" << row << ", value=" << value 
    << (mode == INSERT_VALUES? "(INSERT_VALUES)": "(ADD_VALUES)") << ")";
    
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValue(values_[componentNo], row, value, mode); CHKERRV(ierr);
}

//! set values from another vector
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setValues(PartitionedPetscVec<FunctionSpaceType,nComponents> &rhs)
{
  assert(values_.size() == nComponents);
  assert(!valuesContiguousInUse_);

  VLOG(3) << "\"" << this->name_ << "\" setValues(rhs \"" << rhs.name() << "\")";
  
  // copy the global vectors
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecCopy(rhs.valuesGlobal(componentNo), values_[componentNo]); CHKERRV(ierr);
  }
}

//! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  assert(componentNo >= 0);
  assert(componentNo < values_.size());
  assert(values_.size() == nComponents);
  //assert(!valuesContiguousInUse_); // allow getValues even if valuesContiguousInUse is set

  // this wraps the standard PETSc VecGetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecGetValues(values_[componentNo], ni, ix, y); CHKERRV(ierr); 
  
  // debugging output
  if (VLOG_IS_ON(3))
  {
    std::stringstream str;
    str << "\"" << this->name_ << "\" getValues(componentNo=" << componentNo << ", indices=";
    for (int i = 0; i < ni; i++)
    {
      str << ix[i] << " ";
    }
    str << "): ";
    for (int i = 0; i < ni; i++)
    {
      str << y[i] << " ";
    }
    VLOG(3) << str.str();
  }
}

//! wrapper to the PETSc VecGetValues, on the global vector with global indexing
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
getValuesGlobalPetscIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  getValues(componentNo, ni, ix, y);
}

//! set all entries to zero, wraps VecZeroEntries
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
zeroEntries()
{
  assert(values_.size() == nComponents);
  assert(!valuesContiguousInUse_);

  VLOG(3) << "\"" << this->name_ << "\" zeroEntries";
  
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecZeroEntries(values_[componentNo]); CHKERRV(ierr);
  }
}

//! get the local Vector of a specified component
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
Vec &PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
valuesLocal(int componentNo)
{
  assert(componentNo >= 0);
  assert(componentNo < values_.size());
  assert(values_.size() == nComponents);
  
  return values_[componentNo];
}

//! get the global Vector of a specified component
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
Vec &PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
valuesGlobal(int componentNo)
{
  assert(componentNo >= 0);
  assert(componentNo < values_.size());
  assert(values_.size() == nComponents);
  
  return values_[componentNo];
}

//! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
//! after manipulation of the vector has finished one has to call restoreContiguousValuesGlobal
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
Vec &PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
getContiguousValuesGlobal()
{
  if (nComponents == 1)
  {
    return values_[0];
  }

  PetscErrorCode ierr;

  // create contiguos vector if it does not exist yet
  if (valuesContiguous_ == PETSC_NULL)
  {
    ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &valuesContiguous_); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
    ierr = PetscObjectSetName((PetscObject) valuesContiguous_, this->name_.c_str()); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    // initialize size of vector
    int nEntriesLocal = this->meshPartition_->nDofs() * nComponents;
    int nEntriesGlobal = nEntriesLocal;
    ierr = VecSetSizes(valuesContiguous_, nEntriesLocal, nEntriesGlobal); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(valuesContiguous_); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    LOG(DEBUG) << "\"" << this->name_ << "\" (unstructured) create valuesContiguous_, nComponents = " << nComponents
      << ", nEntriesLocal = " << nEntriesLocal << ", nEntriesGlobal = " << nEntriesGlobal;
  }

  double *valuesDataContiguous;
  ierr = VecGetArray(valuesContiguous_, &valuesDataContiguous); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    const double *valuesDataComponent;
    ierr = VecGetArrayRead(values_[componentNo], &valuesDataComponent); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    VLOG(1) << "  copy " << this->meshPartition_->nDofs()*sizeof(double) << " bytes to contiguous array";
    memcpy(
      valuesDataContiguous + componentNo*this->meshPartition_->nDofs(),
      valuesDataComponent,
      this->meshPartition_->nDofs()*sizeof(double)
    );

    ierr = VecRestoreArrayRead(values_[componentNo], &valuesDataComponent); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
  }

  ierr = VecRestoreArray(valuesContiguous_, &valuesDataContiguous); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
  valuesContiguousInUse_ = true;
  return valuesContiguous_;
}

//! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
//! this has to be called
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
restoreContiguousValuesGlobal()
{
  if (nComponents == 1)
    return;

  assert(valuesContiguous_ != PETSC_NULL);
  if (!valuesContiguousInUse_)
  {
    LOG(FATAL) << "Called restoreContiguousValuesGlobal() without previous getContiguousValuesGlobal()!";
  }

  PetscErrorCode ierr;
  const double *valuesDataContiguous;
  ierr = VecGetArrayRead(valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    double *valuesDataComponent;
    ierr = VecGetArray(values_[componentNo], &valuesDataComponent); CHKERRV(ierr);

    VLOG(1) << "  component " << componentNo << ", copy " << this->meshPartition_->nDofs() << " values, " << this->meshPartition_->nDofs()*sizeof(double) << " bytes from contiguous array";
    memcpy(
      valuesDataComponent,
      valuesDataContiguous + componentNo*this->meshPartition_->nDofs(),
      this->meshPartition_->nDofs()*sizeof(double)
    );

    ierr = VecRestoreArray(values_[componentNo], &valuesDataComponent); CHKERRV(ierr);
  }

  ierr = VecRestoreArrayRead(valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);
  valuesContiguousInUse_ = false;
}

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
output(std::ostream &stream)
{
#ifndef NDEBUG  
  // this method gets all values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo, nRanks;
  MPI_Comm_rank(this->meshPartition_->mpiCommunicator(), &ownRankNo);
  MPI_Comm_size(this->meshPartition_->mpiCommunicator(), &nRanks);
    
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    Vec vector = values_[componentNo];
    if (valuesContiguous_)
      vector = valuesContiguous_;

    // retrieve local values
    int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
    std::vector<double> localValues(nDofsLocal);
    
    VecGetValuesLocal(vector, nDofsLocal, this->meshPartition_->dofNosLocal(), localValues.data());
    
    // gather the local sizes of the vectors to rank 0
    std::vector<int> localSizes(nRanks);
    localSizes[ownRankNo] = nDofsLocal;
    // MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
    // Note that the recvcount argument at the root indicates the number of items it receives from each process, not the total number of items it receives.
    MPI_Gather(localSizes.data() + ownRankNo, 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, this->meshPartition_->mpiCommunicator());
    
    // determine the maximum size/number of vector entries on any rank
    int maxLocalSize;
    MPI_Allreduce(localSizes.data() + ownRankNo, &maxLocalSize, 1, MPI_INT, MPI_MAX, this->meshPartition_->mpiCommunicator());
    
    // gather all values to rank 0
    std::vector<double> recvBuffer(maxLocalSize*nRanks);
    std::vector<double> sendBuffer(maxLocalSize,0.0);
    std::copy(localSizes.begin(), localSizes.end(), sendBuffer.begin());
    
    MPI_Gather(sendBuffer.data(), maxLocalSize, MPI_DOUBLE, recvBuffer.data(), maxLocalSize, MPI_DOUBLE, 0, this->meshPartition_->mpiCommunicator());
    
    if (ownRankNo == 0)
    {
      stream << "vector \"" << this->name_ << "\"" << std::endl;
      stream << "component " << componentNo << ": ";
      for (int rankNo = 0; rankNo < nRanks; rankNo++)
      {
        for (int i = 0; i < localSizes[rankNo]; i++)
        {
          stream << recvBuffer[rankNo*maxLocalSize + i] << "  ";
        }
      }
      stream << std::endl; 
    }

    PetscViewer viewer;
    static int counter = 0;
    std::stringstream vectorOutputFilename;
    vectorOutputFilename << "vector_" << counter++ << ".txt";
    PetscErrorCode ierr;
    ierr = PetscViewerASCIIOpen(this->meshPartition_->mpiCommunicator(), vectorOutputFilename.str().c_str(), &viewer); CHKERRV(ierr);
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INDEX); CHKERRV(ierr);
    ierr = VecView(vector, viewer); CHKERRV(ierr);

    if (ownRankNo == 0)
    {
      stream << "(Vector also written to \"" << vectorOutputFilename.str() << "\".)";
    }
  }
#endif
}

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, PartitionedPetscVec<FunctionSpaceType,nComponents> &vector)
{
  vector.output(stream);
  return stream;
}
