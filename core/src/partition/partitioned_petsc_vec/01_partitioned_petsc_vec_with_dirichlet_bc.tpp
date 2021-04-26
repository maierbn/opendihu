#include "partition/partitioned_petsc_vec/01_partitioned_petsc_vec_with_dirichlet_bc.h"

#include "utility/mpi_utility.h"

//! constructor
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
PartitionedPetscVecWithDirichletBc(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition,
                                   std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,nComponentsDirichletBc>> dirichletBoundaryConditions, std::string name) :
  PartitionedPetscVec<FunctionSpaceType,nComponents>(meshPartition, name), dirichletBoundaryConditions_(dirichletBoundaryConditions)
{
  initialize();
  createVector();
}

//! constructor
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
PartitionedPetscVecWithDirichletBc(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition,
                                   std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,nComponentsDirichletBc>> dirichletBoundaryConditions, std::string name,
                                   bool dummy) :
  PartitionedPetscVec<FunctionSpaceType,nComponents>(meshPartition, name), dirichletBoundaryConditions_(dirichletBoundaryConditions)
{
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
global_no_t PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nonBCDofNoGlobal(int componentNo, dof_no_t localDofNo) const
{
  return dofNoLocalToDofNoNonBcGlobal_[componentNo][localDofNo];
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
dof_no_t PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nonBCDofNoLocal(int componentNo, dof_no_t localDofNo) const
{
  VLOG(2) << "dofNoLocalToDofNoNonBcLocal_[component " << componentNo << "][dofNolocal " << localDofNo << "]=" << dofNoLocalToDofNoNonBcLocal_[componentNo][localDofNo];
  if (dofNoLocalToDofNoNonBcLocal_[componentNo][localDofNo] == -1)
  {
    VLOG(2) << "(=-1)";
  }
  VLOG(2) << "isPrescribed: " << isPrescribed(componentNo, localDofNo);

  return dofNoLocalToDofNoNonBcLocal_[componentNo][localDofNo];
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
bool PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
isPrescribed(int componentNo, dof_no_t localDofNo) const
{
  return isPrescribed_[componentNo][localDofNo];
}

//! get a reference to the internal dofNoLocalToDofNoNonBcGlobal_ data structure
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
const std::array<std::vector<dof_no_t>,nComponents> &PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
dofNoLocalToDofNoNonBcGlobal() const
{
  return dofNoLocalToDofNoNonBcGlobal_;
}

//! get number of local entries
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
int PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nEntriesLocal() const
{
  return nEntriesLocal_;
}

//! get number of global entries
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
global_no_t PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nEntriesGlobal() const
{
  return nEntriesGlobal_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
Vec &PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
valuesGlobal()
{
  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedGlobal)
  {
    VLOG(1) << "valuesGlobal called in not combined-global vector representation ("
      << this->getCurrentRepresentationString()
      <<"), now set to combined-global (without considering ghost dofs, call finishGhostManipulation if ghosts are needed!)";
    setRepresentationGlobal();
  }

  VLOG(2) << "\"" << this->name_ << "\" valuesGlobal()";

  //LOG(DEBUG) << "valuesGlobal, return vectorCombinedWithoutDirichletDofsGlobal_";

  return vectorCombinedWithoutDirichletDofsGlobal_;
}

//! Communicates the ghost values from the global vectors to the local vector and sets the representation to local.
//! The representation has to be global, afterwards it is set to local.
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
startGhostManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" startGhostManipulation";

  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedGlobal)
  {
    LOG(FATAL) << "\"" << this->name_ << "\", startGhostManipulation called when representation is not combined-global (but "
      << this->getCurrentRepresentationString()
      << "), this overwrites the previous values and fetches the last from the global vectors!" << std::endl
      << "Call setRepresentationGlobal() before startGhostManipulation() or check if startGhostManipulation() "
      << "is even necessary (because the representation is already local).";
  }

  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;

  ierr = VecGhostUpdateBegin(vectorCombinedWithoutDirichletDofsGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
  ierr = VecGhostUpdateEnd(vectorCombinedWithoutDirichletDofsGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
  ierr = VecGhostGetLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

  this->currentRepresentation_ = Partition::values_representation_t::representationCombinedLocal;
}

//! Communicates the ghost values from the local vectors back to the global vector and sets the representation to global.
//! The representation has to be local, afterwards it is set to global.
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
finishGhostManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" finishGhostManipulation";

  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedLocal)
  {
    LOG(ERROR) << "\"" << this->name_ << "\", finishGhostManipulation called when representation is not combined-local (it is "
      << this->getCurrentRepresentationString()
      << "), (probably no previous startGhostManipulation)";
  }

  // Copy the local values vectors into the global vector. ADD_VALUES means that ghost values are reduced (summed up)
  PetscErrorCode ierr;
  ierr = VecGhostRestoreLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);
  ierr = VecGhostUpdateBegin(vectorCombinedWithoutDirichletDofsGlobal_, ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);
  ierr = VecGhostUpdateEnd(vectorCombinedWithoutDirichletDofsGlobal_, ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);

  this->currentRepresentation_ = Partition::values_representation_t::representationCombinedGlobal;
}

// set the internal representation to be global, i.e. using the global vectors
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setRepresentationGlobal()
{
  VLOG(2) << "\"" << this->name_ << "\" PartitionedPetscVecWithDirichletBc::setRepresentationGlobal, previous representation: "
    << this->getCurrentRepresentationString();

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    return;
  }

    // set the internal representation to be global, i.e. using the global vectors
  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    // restore local vector back into global vector, does not consider ghost dofs
    // use finishGhostManipulation to handle ghost dofs, that requires communication

    PetscErrorCode ierr;
    ierr = VecGhostRestoreLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

    this->currentRepresentation_ = Partition::values_representation_t::representationCombinedGlobal;
  }
  else
  {
    LOG(FATAL) << "Cannot set vector representation from \"" << this->getCurrentRepresentationString()
      << "\" to \"combined-global\".";
  }
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setRepresentationLocal()
{
  VLOG(2) << "\"" << this->name_ << "\" setRepresentationLocal, previous representation: "
    << this->getCurrentRepresentationString();

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    // retrieve the local vector from the global vector, does not consider ghost dofs
    // use startGhostManipulation to also set the correct values for the ghost dofs, that requires communication

    PetscErrorCode ierr;
    ierr = VecGhostGetLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

    this->currentRepresentation_ = Partition::values_representation_t::representationCombinedLocal;
  }
  else
  {
    LOG(FATAL) << "Cannot set vector representation from \"" << this->getCurrentRepresentationString()
      << "\" to \"combined-local\".";
  }
}

//! zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
zeroGhostBuffer()
{
  VLOG(2) << "\"" << this->name_ << "\" zeroGhostBuffer";

  // set local ghost values to 0
  std::vector<dof_no_t> indices(nNonBcDofsGhosts_);
  std::iota(indices.begin(), indices.end(), nEntriesLocal_);
  std::vector<double> values(nNonBcDofsGhosts_, 0.0);

  PetscErrorCode ierr;
  ierr = VecSetValues(vectorCombinedWithoutDirichletDofsLocal_, nNonBcDofsGhosts_, indices.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
}

//! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  VLOG(1) << "\"" << this->name_ << "\" setValues, representation: " << this->getCurrentRepresentationString();

  assert(componentNo >= 0 && componentNo < nComponents);

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    // determine new indices
    std::vector<dof_no_t> indices;
    std::vector<double> values;

    indices.reserve(ni);
    values.reserve(ni);

  PetscErrorCode ierr;
#ifndef NDEBUG
    dof_no_t ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vectorCombinedWithoutDirichletDofsGlobal_, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
#endif

    for (dof_no_t i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithoutGhosts());

      if (!isPrescribed_[componentNo][ix[i]])
      {
        dof_no_t nonBcIndexGlobal = nonBCDofNoGlobal(componentNo, ix[i]);
#ifndef NDEBUG
        if (!(nonBcIndexGlobal >= ownershipBegin && nonBcIndexGlobal < ownershipEnd))
        {
          LOG(ERROR) << "setValues \"" << this->name_ << "\", i=" << i << ", ix[i]=" << ix[i] << ", nonBcIndexGlobal=" << nonBcIndexGlobal
            << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
        }
        assert(nonBcIndexGlobal >= ownershipBegin && nonBcIndexGlobal < ownershipEnd);
#endif

        indices.push_back(nonBcIndexGlobal);
        values.push_back(y[i]);
      }
    }

    VLOG(2) << "non-bc indices: " << indices << ", values: " << values;
    // this wraps the standard PETSc VecSetValues on the local vector
    ierr = VecSetValues(vectorCombinedWithoutDirichletDofsGlobal_, indices.size(), indices.data(), values.data(), iora); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
  else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    // determine new indices
    std::vector<dof_no_t> indices;
    std::vector<double> values;

    indices.reserve(ni);
    values.reserve(ni);

    for (dof_no_t i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithGhosts());

      if (!isPrescribed_[componentNo][ix[i]])
      {
        indices.push_back(nonBCDofNoLocal(componentNo, ix[i]));
        values.push_back(y[i]);
      }
    }

    VLOG(2) << "non-bc indices: " << indices << ", values: " << values;
    // this wraps the standard PETSc VecSetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValues(vectorCombinedWithoutDirichletDofsLocal_, indices.size(), indices.data(), values.data(), iora); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
}

//! wrapper to the PETSc VecSetValue, acting only on the local data
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  VLOG(1) << "\"" << this->name_ << "\" setValue row=" << row << ", representation: " << this->getCurrentRepresentationString();

  assert(componentNo >= 0 && componentNo < nComponents);

  // replace dirichlet BC values with the prescribed values
  if (isPrescribed_[componentNo][row])
  {
    VLOG(2) << "row " << row << ", value is prescribed, do not change";
    return;
  }

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    assert(row < this->meshPartition_->nDofsLocalWithoutGhosts());

  PetscErrorCode ierr;
#ifndef NDEBUG
    dof_no_t ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vectorCombinedWithoutDirichletDofsGlobal_, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
#endif

    // determine new index
    row = nonBCDofNoGlobal(componentNo, row);

#ifndef NDEBUG
    if (!(row >= ownershipBegin && row < ownershipEnd))
    {
      LOG(ERROR) << "setValue \"" << this->name_ << "\", row=" << row
        << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
    }
    assert(row >= ownershipBegin && row < ownershipEnd);
#endif

    VLOG(2) << "set value " << value << " at non-bc global " << row;

    // this wraps the standard PETSc VecSetValues on the local vector
    ierr = VecSetValue(vectorCombinedWithoutDirichletDofsGlobal_, row, value, mode); CHKERRV(ierr);
  }
  else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    assert(row < this->meshPartition_->nDofsLocalWithGhosts());

    VLOG(2) << "set value at dof local " << row << " (n without ghosts: "
      << this->meshPartition_->nDofsLocalWithoutGhosts() << ", with: " << this->meshPartition_->nDofsLocalWithGhosts() << ")";

    // determine new index
    row = nonBCDofNoLocal(componentNo, row);

    VLOG(2) << "set value " << value << " at non-bc local " << row;

    assert(row < nEntriesLocal_ + nNonBcDofsGhosts_);

    // this wraps the standard PETSc VecSetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValue(vectorCombinedWithoutDirichletDofsLocal_, row, value, mode); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
}

//! wrapper to the PETSc VecSetValue, acting on the local data or global data, the row is local dof no
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setValue(int componentNo, Vc::int_v row, Vc::double_v value, InsertMode mode)
{
  // loop over rows and set non-vectorized values
  for (int rowIndex = 0; rowIndex < Vc::double_v::size(); rowIndex++)
  {
    if (row[rowIndex] == -1)
      break;
    this->setValue(componentNo, (int)(row[rowIndex]), (double)(value[rowIndex]), mode);
  }
}

//! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const
{
  VLOG(1) << "\"" << this->name_ << "\" getValues, representation: " << this->getCurrentRepresentationString();

  if (componentNo >= nComponents)
    LOG(FATAL) << "\"" << this->name_ << "\" getValues, representation: " << this->getCurrentRepresentationString()
      << ", componentNo is invalid, componentNo " << componentNo << ", nComponents: " << nComponents;

  assert(componentNo >= 0 && componentNo < nComponents);

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {

    PetscErrorCode ierr;
#ifndef NDEBUG
    dof_no_t ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vectorCombinedWithoutDirichletDofsGlobal_, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);

    // check that indices are in range
    for (dof_no_t i = 0; i < ni; i++)
    {
      if (ix[i] >= this->meshPartition_->nDofsLocalWithoutGhosts())
      {
        LOG(ERROR) << "ni = " << ni << ", ix[" << i << "]=" << ix[i] << " >= " << this->meshPartition_->nDofsLocalWithoutGhosts();
      }
    }
#endif

    // determine new global indices
    std::vector<dof_no_t> indices(ni);
    for (dof_no_t i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithoutGhosts());
      indices[i] = nonBCDofNoGlobal(componentNo, ix[i]);

      // check if indices are in range that is owned by own rank
#ifndef NDEBUG
      if (!(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd)))
      {
        LOG(ERROR) << "getValues \"" << this->name_ << "\", i=" << i << ", ix[i]=" << ix[i] << ", indices[i]=" << indices[i]
          << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
      }
      assert(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd));
#endif
    }

    // this wraps the standard PETSc VecGetValues on the local vector
    ierr = VecGetValues(vectorCombinedWithoutDirichletDofsGlobal_, ni, indices.data(), y); CHKERRV(ierr);

#ifndef NDEBUG
    VLOG(1) << "non-bc global indices: " << indices;
    VLOG(1) << "got values: ";
    for (dof_no_t i = 0; i < ni; i++)
      VLOG(1) << y[i];
#endif

    // replace dirichlet BC values with the prescribed values
    for (dof_no_t i = 0; i < ni; i++)
    {
      if (isPrescribed_[componentNo][ix[i]] && componentNo < nComponentsDirichletBc)
      {
        y[i] = boundaryConditionValues_[componentNo][ix[i]];
      }
    }

#ifndef NDEBUG
    VLOG(1) << "modified values: ";
    for (dof_no_t i = 0; i < ni; i++)
      VLOG(1) << y[i];

    VLOG(1) << "component " << componentNo << ", nComponentsDirichletBc: " << nComponentsDirichletBc
      << "\n isPrescribed: " << isPrescribed_[componentNo]
      << "\n boundaryConditionValues: " << boundaryConditionValues_[componentNo] << "\nmodified values: ";
    for (dof_no_t i = 0; i < ni; i++)
      VLOG(1) << "local index " << ix[i] << " global " << indices[i] << ", isPrescribed: " << isPrescribed_[componentNo][ix[i]] << ", bc value: " << boundaryConditionValues_[componentNo][ix[i]] << ", final value: " << y[i];
#endif
  }
  else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    // determine new indices
    std::vector<dof_no_t> indices(ni);
    for (dof_no_t i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithGhosts());
      indices[i] = nonBCDofNoLocal(componentNo, ix[i]);
    }

    // this wraps the standard PETSc VecGetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecGetValues(vectorCombinedWithoutDirichletDofsLocal_, ni, indices.data(), y); CHKERRV(ierr);

#ifndef NDEBUG
    VLOG(1) << "non-bc local indices: " << indices;
    VLOG(1) << "got values: ";
    for (dof_no_t i = 0; i < ni; i++)
      VLOG(1) << y[i];
#endif

    // replace dirichlet BC values with the prescribed values
    for (dof_no_t i = 0; i < ni; i++)
    {
      if (isPrescribed_[componentNo][ix[i]] && componentNo < nComponentsDirichletBc)
      {
        y[i] = boundaryConditionValues_[componentNo][ix[i]];
      }
    }

#ifndef NDEBUG
    VLOG(1) << "component " << componentNo << ", nComponentsDirichletBc: " << nComponentsDirichletBc
      << "\n isPrescribed: " << isPrescribed_[componentNo]
      << "\n boundaryConditionValues: " << boundaryConditionValues_[componentNo] << "\nmodified values: ";
    for (dof_no_t i = 0; i < ni; i++)
      VLOG(1) << "local index " << ix[i] << ", value: " << y[i];
#endif
  }
}

//! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
getValuesGlobal(Vec vector, int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const
{
  LOG(FATAL) << "getValuesGlobal should not be called.";

  PetscErrorCode ierr;
  std::string name;
  char *cName;
  ierr = PetscObjectGetName((PetscObject)vector, (const char **)&cName); CHKERRV(ierr);
  name = cName;

  VLOG(3) << " getValuesGlobal(\"" << name << "\")";

#ifndef NDEBUG
    dof_no_t ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vector, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
#endif

  // determine new global indices
  std::vector<dof_no_t> indices(ni);
  for (dof_no_t i = 0; i < ni; i++)
  {
    assert(ix[i] < this->meshPartition_->nDofsLocalWithoutGhosts());
    indices[i] = nonBCDofNoGlobal(componentNo, ix[i]);

#ifndef NDEBUG
    if (!(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd)))
    {
      LOG(ERROR) << "getValuesGlobal \"" << this->name_ << "\", i=" << i << ", ix[i]=" << ix[i] << ", indices[i]=" << indices[i]
        << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
    }
    assert(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd));
#endif
  }

  VLOG(1) << "non-bc global indices: " << indices;

  // this wraps the standard PETSc VecGetValues on the local vector
  ierr = VecGetValues(vector, ni, indices.data(), y); CHKERRV(ierr);

  VLOG(1) << "got values: ";
  for (dof_no_t i = 0; i < ni; i++)
    VLOG(1) << y[i];

  // replace dirichlet BC values with the prescribed values
  for (dof_no_t i = 0; i < ni; i++)
  {
    if (isPrescribed_[componentNo][ix[i]])
    {
      if (componentNo < nComponentsDirichletBc)
      {
        y[i] = boundaryConditionValues_[componentNo][ix[i]];
      }
    }
  }

  VLOG(1) << "modified values: ";
  for (dof_no_t i = 0; i < ni; i++)
    VLOG(1) << y[i];
}

//! get a single value
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
double PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
getValue(int componentNo, PetscInt row) const
{
  double result;
  this->getValues(componentNo, 1, &row, &result);
  return result;
}

//! set all entries to zero, wraps VecZeroEntries
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
zeroEntries()
{
  VLOG(3) << "\"" << this->name_ << "\" zeroEntries, representation: " << this->getCurrentRepresentationString();

  PetscErrorCode ierr;
  //if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    ierr = VecZeroEntries(vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);
  }
  //else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    ierr = VecZeroEntries(vectorCombinedWithoutDirichletDofsGlobal_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
updateDirichletBoundaryConditions(const std::vector<std::pair<global_no_t,std::array<double,nComponentsDirichletBc>>> &newValues,
                                  bool inputMeshIsGlobal)
{
  LOG(DEBUG) << "updateDirichletBoundaryConditions, newValues: " << newValues;

  for (int itemIndex = 0; itemIndex < newValues.size(); itemIndex++)
  {
    dof_no_t dofNoLocal = newValues[itemIndex].first;

    // transform global dof no to local dof no if needed
    if (inputMeshIsGlobal)
    {
      global_no_t dofNoGlobal = newValues[itemIndex].first;

      bool isLocal = false;
      dofNoLocal = this->meshPartition_->getDofNoLocal(dofNoGlobal, isLocal);

      if (!isLocal)
        continue;
    }

    // assign new values
    std::array<double,nComponentsDirichletBc> values = newValues[itemIndex].second;

    for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
    {
      boundaryConditionValues_[componentNo][dofNoLocal] = values[componentNo];
    }
  }

  communicateBoundaryConditionGhostValues();
}


template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
bool PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
containsNanOrInf()
{
  // get all local values
  std::vector<double> values;

  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    values.resize(this->meshPartition_->nDofsLocalWithoutGhosts());

    // void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const;
    this->getValues(componentNo, this->meshPartition_->nDofsLocalWithoutGhosts(), this->meshPartition_->dofNosLocal().data(), values.data());

    // loop over values and check if they are neither nan nor inf
    for (int i = 0; i < values.size(); i++)
    {
      if (!std::isfinite(values[i]) || fabs(values[i]) > 1e+75)
      {
        LOG(ERROR) << "containsNanOrInf(): value " << i << "/" << values.size() << ", component " << componentNo << "/" << nComponents << ": " << values[i];
        return true;
      }
    }
  }
  return false;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
Vec &PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
valuesGlobalReference()
{
  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedGlobal)
  {
    LOG(FATAL) << "Calling valuesGlobalReference() when representation is not combined-global, but " << this->getCurrentRepresentationString();
  }
  return vectorCombinedWithoutDirichletDofsGlobal_;
}

//! output the vector to stream, for debugging
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
output(std::ostream &stream) const
{
  stream << "\"" << this->name_ << "\" " << this->nEntriesGlobal_ << " entries global, " << nEntriesLocal_ << " entries local." << std::endl
    << "nNonBcDofsWithoutGhosts: " << nNonBcDofsWithoutGhosts_ << ", nNonBcDofsGhosts: " << nNonBcDofsGhosts_ << ", nonBcGhostDofNosGlobal_: " << nonBcGhostDofNosGlobal_ << std::endl
    << "dofNoLocalToDofNoNonBcGlobal_: " << dofNoLocalToDofNoNonBcGlobal_ << std::endl
    << "dofNoLocalToDofNoNonBcLocal_:  " << dofNoLocalToDofNoNonBcLocal_ << std::endl
    << "boundaryConditionValues_: " << boundaryConditionValues_ << ", isPrescribed_: " << isPrescribed_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscVecWithDirichletBc<FunctionSpaceType,nComponents> &vector)
{
  vector.output(stream);
  return stream;
}
