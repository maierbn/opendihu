#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
PartitionedPetscVec<FunctionSpaceType,nComponents,DummyForTraits>::
PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, std::string name) :
  PartitionedPetscVecBase<FunctionSpaceType>(meshPartition, name)
{
  createVector();
}

//! constructor from existing vector, values are not copied
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
template<int nComponents2>
PartitionedPetscVec<FunctionSpaceType,nComponents,DummyForTraits>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpaceType,nComponents2> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpaceType>(rhs.meshPartition(), name)
 
{
  createVector();

  setValues(rhs);
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

// set the internal representation to be global, i.e. using the global vectors
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setRepresentationGlobal()
{
  if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
  {
    // this sets the representation to local
    restoreValuesContiguous();
  }
  this->currentRepresentation_ = Partition::values_representation_t::representationGlobal;
}

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setRepresentationLocal()
{
  // global and local are the same for unstructured meshes
  setRepresentationGlobal();

  this->currentRepresentation_ = Partition::values_representation_t::representationLocal;
}

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setRepresentationContiguous()
{
  if (this->currentRepresentation_ == Partition::values_representation_t::representationLocal)
  {
    getValuesContiguous();
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
    str << "\"" << this->name_ << "\", representation \""
      << Partition::valuesRepresentationString[this->currentRepresentation_]
      << "\", setValues(componentNo=" << componentNo << ", indices=";
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
  
  if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
  {
    if (componentNo > 0)
    {
      // shift indices
      std::vector<PetscInt> indices(ni);
      for (int i = 0; i < ni; i++)
      {
        indices[i] = ix[i] + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts();
      }

      // this wraps the standard PETSc VecSetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecSetValues(valuesContiguous_, ni, indices.data(), y, iora); CHKERRV(ierr);
    }
    else
    {
      // this wraps the standard PETSc VecSetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecSetValues(valuesContiguous_, ni, ix, y, iora); CHKERRV(ierr);
    }
  }
  else
  {
    PetscErrorCode ierr;
    ierr = VecSetValues(values_[componentNo], ni, ix, y, iora); CHKERRV(ierr);

    // get value
    double value;
    ierr = VecGetValues(values_[componentNo], ni, ix, &value); CHKERRV(ierr);
    LOG(DEBUG) << "retrieved value: " << value;
  }
}

//! wrapper to the PETSc VecSetValue, acting only on the local data
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  assert(componentNo >= 0);
  assert(componentNo < values_.size());
  assert(values_.size() == nComponents);

  // debugging output
  VLOG(3) << "\"" << this->name_ << "\" setValue(componentNo=" << componentNo << ", row=" << row << ", value=" << value 
    << (mode == INSERT_VALUES? "(INSERT_VALUES)": "(ADD_VALUES)") << ")";
    
  if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
  {
    PetscErrorCode ierr;
    PetscInt index = row + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts();
    ierr = VecSetValue(valuesContiguous_, index, value, mode); CHKERRV(ierr);
  }
  else
  {
    // this wraps the standard PETSc VecSetValue on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValue(values_[componentNo], row, value, mode); CHKERRV(ierr);
  }
}

//! set values from another vector
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
template<int nComponents2>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
setValues(PartitionedPetscVec<FunctionSpaceType,nComponents2> &rhs)
{
  assert(values_.size() == nComponents);

  VLOG(3) << "\"" << this->name_ << "\" setValues(rhs \"" << rhs.name() << "\"), rhs representation: "
    << Partition::valuesRepresentationString[rhs.currentRepresentation()];
  
  assert(this->currentRepresentation_ != Partition::values_representation_t::representationContiguous);

  // copy existing values
  PetscErrorCode ierr;
  if (rhs.currentRepresentation() == Partition::values_representation_t::representationContiguous)
  {
    ierr = VecCopy(rhs.getValuesContiguous(), valuesContiguous_); CHKERRV(ierr);
  }
  else
  {
    for (int componentNo = 0; componentNo < std::min(nComponents,nComponents2); componentNo++)
    {
      ierr = VecCopy(rhs.valuesGlobal(componentNo), values_[componentNo]); CHKERRV(ierr);
    }
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

  // depending on which data representation is active, use vectorLocal or valuesContiguous
  if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
  {
    if (componentNo > 0)
    {
      // shift indices
      std::vector<PetscInt> indices(ni);
      for (int i = 0; i < ni; i++)
      {
        indices[i] = ix[i] + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts();
      }

      // this wraps the standard PETSc VecGetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecGetValues(valuesContiguous_, ni, indices.data(), y); CHKERRV(ierr);
    }
    else
    {
      // this wraps the standard PETSc VecGetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecGetValues(valuesContiguous_, ni, ix, y); CHKERRV(ierr);
    }
  }
  else
  {
    // this wraps the standard PETSc VecGetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecGetValues(values_[componentNo], ni, ix, y); CHKERRV(ierr);
  }

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

  VLOG(3) << "\"" << this->name_ << "\" zeroEntries";
  
  PetscErrorCode ierr;
  if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
  {
    ierr = VecZeroEntries(valuesContiguous_); CHKERRV(ierr);
  }
  else
  {
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      ierr = VecZeroEntries(values_[componentNo]); CHKERRV(ierr);
    }
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
//! after manipulation of the vector has finished one has to call restoreValuesContiguous
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
Vec &PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
getValuesContiguous()
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
  this->currentRepresentation_ = Partition::values_representation_t::representationContiguous;
  return valuesContiguous_;
}

//! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
//! this has to be called
template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
restoreValuesContiguous()
{
  if (nComponents == 1)
    return;

  assert(valuesContiguous_ != PETSC_NULL);
  if (this->currentRepresentation_ != Partition::values_representation_t::representationContiguous)
  {
    LOG(FATAL) << "Called restoreValuesContiguous() without previous getValuesContiguous()!";
  }

  PetscErrorCode ierr;
  const double *valuesDataContiguous;
  ierr = VecGetArrayRead(valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    double *valuesDataComponent;
    ierr = VecGetArray(values_[componentNo], &valuesDataComponent); CHKERRV(ierr);

    VLOG(1) << "  \"" << this->name_ << "\", component " << componentNo << ", copy " << this->meshPartition_->nDofs() << " values, " << this->meshPartition_->nDofs()*sizeof(double) << " bytes from contiguous array";
    memcpy(
      valuesDataComponent,
      valuesDataContiguous + componentNo*this->meshPartition_->nDofs(),
      this->meshPartition_->nDofs()*sizeof(double)
    );

    ierr = VecRestoreArray(values_[componentNo], &valuesDataComponent); CHKERRV(ierr);
  }

  ierr = VecRestoreArrayRead(valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);
  this->currentRepresentation_ = Partition::values_representation_t::representationLocal;
}

template<typename FunctionSpaceType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<FunctionSpaceType, nComponents, DummyForTraits>::
output(std::ostream &stream)
{
#ifndef NDEBUG  
  // this method gets all values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo = this->meshPartition_->ownRankNo();
  PetscMPIInt nRanks = this->meshPartition_->nRanks();

  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    Vec vector = values_[componentNo];
    if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
      vector = valuesContiguous_;

    // retrieve local values
    int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
    std::vector<PetscInt> indices(nDofsLocal);
    if (this->currentRepresentation_ == Partition::values_representation_t::representationContiguous)
    {
      for (int i = 0; i < nDofsLocal; i++)
      {
        indices[i] = this->meshPartition_->dofNosLocal()[i] + componentNo*nDofsLocal;
      }
    }
    else
    {
      for (int i = 0; i < nDofsLocal; i++)
      {
        indices[i] = this->meshPartition_->dofNosLocal()[i];
      }
    }
    std::vector<double> localValues(nDofsLocal);
    PetscErrorCode ierr;
    ierr = VecGetValues(vector, nDofsLocal, indices.data(), localValues.data()); CHKERRV(ierr);
    
    if (ownRankNo == 0)
    {
      stream << "vector \"" << this->name_ << "\", representation \""
        << Partition::valuesRepresentationString[this->currentRepresentation_] << "\"" << std::endl;
      stream << "component " << componentNo << ": [";
      for (int rankNo = 0; rankNo < nRanks; rankNo++)
      {
        for (int i = 0; i < nDofsLocal; i++)
        {
          stream << localValues[i] << "  ";
        }
      }
      stream << "]" << std::endl;
    }

    PetscViewer viewer;
    static int counter = 0;
    std::stringstream vectorOutputFilename;
    vectorOutputFilename << "vector_" << counter++ << ".txt";
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
