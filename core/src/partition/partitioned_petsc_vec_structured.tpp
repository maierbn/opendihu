#include "partition/partitioned_petsc_vec.h"

//! constructor
template<typename MeshType,typename BasisFunctionType,int nComponents>
PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition,
                    std::string name) :
  PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, name), ghostManipulationStarted_(false)
{
  //dm_ = meshPartition->dmElements();
  
  createVector();
}

//! constructor, copy from existing vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(rhs.meshPartition(), name), ghostManipulationStarted_(false)
{
  createVector();

  // copy existing values
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecCopy(rhs.valuesLocal(componentNo), vectorLocal_[componentNo]); CHKERRV(ierr);
    /*ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = VecCopy(rhs.valuesGlobal(componentNo), vectorGlobal_[i]); CHKERRV(ierr);
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);*/
  }
}
  
//! constructor, copy from existing vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
template<int nComponents2>
PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents2> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(rhs.meshPartition(), name), ghostManipulationStarted_(false)
{
  //dm_ = rhs.dm_;
  
  createVector();
  
  // copy existing values
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < std::min(nComponents,nComponents2); componentNo++)
  {
    ierr = VecCopy(rhs.valuesLocal(componentNo), vectorLocal_[componentNo]); CHKERRV(ierr);
    /*ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = VecCopy(rhs.valuesGlobal(componentNo), vectorGlobal_[i]); CHKERRV(ierr);
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);*/
  }
}
  
//! create a distributed Petsc vector, according to partition
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
createVector()
{
  VLOG(2) << "\"" << this->name_ << "\" createVector with " << nComponents << " components, size local: " << this->meshPartition_->nNodesLocalWithoutGhosts() 
    << ", global: " << this->meshPartition_->nNodesGlobal() << ", ghost dof nos global/petsc: " << this->meshPartition_->ghostDofNosGlobalPetsc();
  PetscErrorCode ierr;
  
  // The local vector contains the nodal/dof values for the local portion of the current rank. This includes ghost nodes.
  // The global vector manages the whole data and also only stores the local portion of it on the current rank.
  // It shares the memory with the local vector. The local vector uses local indices whereas the global vector is accessed using globalPetsc indexing.
  // If we want to manipulate data, we have to ensure consistency in the parallel global vector (VecGhostUpdateBegin,VecGhostUpdateEnd)
  // and then fetch the local portion of the global vector together with ghost values into a local vector (VecGhostGetLocalForm),
  // then manipulate the values (using standard VecSetValues/VecGetValues routines) and commit the results (VecGhostRestoreLocalForm, VecGhostUpdateBegin, VecGhostUpdateEnd).
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    const dof_no_t nGhostDofs = this->meshPartition_->nDofsLocalWithGhosts() - this->meshPartition_->nDofsLocalWithoutGhosts();
    ierr = VecCreateGhost(this->meshPartition_->mpiCommunicator(), this->meshPartition_->nDofsLocalWithoutGhosts(),
                          this->meshPartition_->nDofsGlobal(), nGhostDofs, this->meshPartition_->ghostDofNosGlobalPetsc().data(), &vectorGlobal_[componentNo]); CHKERRV(ierr);
    
    // initialize PETSc vector object
    //ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject) vectorGlobal_[componentNo], this->name_.c_str()); CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(vectorGlobal_[componentNo]); CHKERRV(ierr);
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
  }

  // createVector acts like startGhostManipulation as it also gets the local vector (VecGhostGetLocalForm) to work on.
  ghostManipulationStarted_ = true;
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
startGhostManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" startGhostManipulation";
  
  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;
  /*if (ghostManipulationStarted_)
  {
    LOG(DEBUG) << "\"" << this->name_ << "\", startGhostManipulation called multiple times without finishGhostManipulation in between. (Creation of the vector also calls startGhostManipulation). "
      << "That means that locally stored ghost values got overwritten with the actual values from the foreign owner rank.";
  }*/
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecGhostUpdateBegin(vectorGlobal_[componentNo], INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
    //ierr = VecZeroEntries(vectorLocal_[componentNo]); CHKERRV(ierr);
    //ierr = DMGlobalToLocalBegin(*dm_, vectorGlobal_[componentNo], INSERT_VALUES, vectorLocal_[componentNo]); CHKERRV(ierr);
  }
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecGhostUpdateEnd(vectorGlobal_[componentNo], INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
    //ierr = DMGlobalToLocalEnd(*dm_, vectorGlobal_[componentNo], INSERT_VALUES, vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
  }
  ghostManipulationStarted_ = true;
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
finishGhostManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" finishGhostManipulation";
  
  if (!ghostManipulationStarted_)
  {
    LOG(ERROR) << "\"" << this->name_ << "\", finishGhostManipulation called without previous startGhostManipulation";
  }
  
  // Copy the local values vectors into the global vector. ADD_VALUES means that ghost values are reduced (summed up)
  PetscErrorCode ierr;
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
    
    ierr = VecGhostUpdateBegin(vectorGlobal_[componentNo], ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);
    //ierr = VecZeroEntries(vectorGlobal_[componentNo]); CHKERRV(ierr);
    //ierr = DMLocalToGlobalBegin(*dm_, vectorLocal_[componentNo], ADD_VALUES, vectorGlobal_[componentNo]); CHKERRV(ierr);
  }
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecGhostUpdateEnd(vectorGlobal_[componentNo], ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);
    //ierr = DMLocalToGlobalEnd(*dm_, vectorLocal_[componentNo], ADD_VALUES, vectorGlobal_[componentNo]); CHKERRV(ierr);
  }
  ghostManipulationStarted_ = false;
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
zeroGhostBuffer()
{
  VLOG(2) << "\"" << this->name_ << "\" zeroGhostBuffer";

  // set local ghost values to 0
  int nValues = this->meshPartition_->nDofsLocalWithGhosts() - this->meshPartition_->nDofsLocalWithoutGhosts();
  const PetscInt *indices = this->meshPartition_->dofNosLocal().data() + this->meshPartition_->nDofsLocalWithoutGhosts();
  std::vector<double> values(nValues, 0.0);

  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecSetValues(vectorLocal_[componentNo], nValues, indices, values.data(), INSERT_VALUES); CHKERRV(ierr);
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // depending on which data representation is active, use vectorLocal or valuesContiguous
  if (valuesContiguousInUse_)
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
    ierr = VecGetValues(vectorLocal_[componentNo], ni, ix, y); CHKERRV(ierr);
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
    str << ") [valuesContiguousInUse_=" << valuesContiguousInUse_ << "]: ";
    for (int i = 0; i < ni; i++)
    {
      str << y[i] << " ";
    }
    str << "]";
    VLOG(3) << str.str();
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValuesGlobalPetscIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  assert(!valuesContiguousInUse_);

  // this wraps the standard PETSc VecGetValues on the global vector
  PetscErrorCode ierr;
  ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
  ierr = VecGetValues(vectorGlobal_[componentNo], ni, ix, y); CHKERRV(ierr);
  ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
  
  // debugging output
  if (VLOG_IS_ON(3))
  {
    std::stringstream str;
    str << "\"" << this->name_ << "\" getValuesGlobalPetscIndexing(componentNo=" << componentNo << ", indices=";
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

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  if (VLOG_IS_ON(3))
  {
    std::stringstream str;
    str << "\"" << this->name_ << "\" setValues(componentNo=" << componentNo << ", indices=";
    for (int i = 0; i < ni; i++)
    {
      str << ix[i] << " ";
    }
    str << ", values=";
    for (int i = 0; i < ni; i++)
    {
      str << y[i] << " ";
    }
    str << (iora == INSERT_VALUES? "INSERT_VALUES" : (iora == ADD_VALUES? "ADD_VALUES" : "unknown"));
    str << ") " << (valuesContiguousInUse_? " valuesContiguousInUse" : "");
    VLOG(3) << str.str();
  }

  // depending on which data representation is active, use vectorLocal or valuesContiguous
  if (valuesContiguousInUse_)
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

#ifndef NDEBUG

    for (int i = 0; i < ni; i++)
    {
      if (ix[i] < 0 || ix[i] >= this->meshPartition_->nDofsLocalWithGhosts())
      {
        LOG(ERROR) << "\"" << this->name_ << "\" setValues ix[" << i << "]=" << ix[i] << ", nDofsLocalWithGhosts: " << this->meshPartition_->nDofsLocalWithGhosts();
      }
      assert(ix[i] >= 0 && ix[i] < this->meshPartition_->nDofsLocalWithGhosts());
    }
#endif

    // this wraps the standard PETSc VecSetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValues(vectorLocal_[componentNo], ni, ix, y, iora); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  VLOG(3) << "\"" << this->name_ << "\" setValue(componentNo=" << componentNo << ", row=" << row << ", value=" << value
    << (mode == INSERT_VALUES? "INSERT_VALUES" : (mode == ADD_VALUES? "ADD_VALUES" : "unknown"));
  assert(!valuesContiguousInUse_);

  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValue(vectorLocal_[componentNo], row, value, mode); CHKERRV(ierr);
}

//! set values from another vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValues(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents> &rhs)
{
  VLOG(3) << "\"" << this->name_ << "\" setValues(rhs \"" << rhs.name() << "\"), this calls startGhostManipulation()";
  assert(!valuesContiguousInUse_);

  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecCopy(rhs.valuesGlobal(componentNo), vectorGlobal_[componentNo]); CHKERRV(ierr);
  }
  
  startGhostManipulation();
}

/*
//! for a single component vector set all values
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValuesWithGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->meshPartition_->nDofsLocalWithGhosts());
 
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(vectorLocal_[componentNo], values.size(), this->meshPartition_->dofNosLocal().data(), values.data(), petscInsertMode); CHKERRV(ierr);
}

//! for a single component vector set all values, input does not contain values for ghosts
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValuesWithoutGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->meshPartition_->nDofsLocalWithoutGhosts());
 
  // copy all values to a bigger vector that is initialized with 0
  std::vector<double> valuesWithGhosts(this->meshPartition_->nDofsLocalWithGhosts(), 0.0);
  dof_no_t valueIndex = 0;
  for (std::vector<dof_no_t>::const_iterator localDof = this->meshPartition_->nonGhostDofsBegin(); localDof != this->meshPartition_->nonGhostDofsEnd(); localDof++)
  {
    valuesWithGhosts[*localDof] = values[valueIndex++];
  }
  
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(vectorLocal_[componentNo], valuesWithGhosts.size(), this->meshPartition_->localDofNos().data(), valuesWithGhosts.data(), petscInsertMode); CHKERRV(ierr);
}
  */
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
zeroEntries()
{
  VLOG(3) << "\"" << this->name_ << "\" zeroEntries";
  
  assert(!valuesContiguousInUse_);

  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecZeroEntries(vectorLocal_[componentNo]); CHKERRV(ierr);
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
Vec &PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
valuesLocal(int componentNo)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  return vectorLocal_[componentNo];
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
Vec &PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
valuesGlobal(int componentNo)
{
  assert(componentNo >= 0 && componentNo < nComponents);

  // We still may need to call VecGhostRestoreLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]) here?
  // Then after operation on the global vector is done, do we need to call VecGhostGetLocalForm again?
  return vectorGlobal_[componentNo];
}

//! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
//! after manipulation of the vector has finished one has to call restoreContiguousValuesGlobal
template<typename MeshType,typename BasisFunctionType,int nComponents>
Vec &PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getContiguousValuesGlobal()
{
  if (nComponents == 1)
  {
    return vectorGlobal_[0];
  }

  // if the contiguous representation is already being used, return contiguous vector
  if (valuesContiguousInUse_)
  {
    return valuesContiguous_;
  }

  PetscErrorCode ierr;

  // create contiguos vector if it does not exist yet
  if (valuesContiguous_ == PETSC_NULL)
  {
    ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &valuesContiguous_); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
    ierr = PetscObjectSetName((PetscObject) valuesContiguous_, this->name_.c_str()); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    // initialize size of vector
    int nEntriesLocal = this->meshPartition_->nDofsLocalWithoutGhosts() * nComponents;
    int nEntriesGlobal = nEntriesLocal;
    ierr = VecSetSizes(valuesContiguous_, nEntriesLocal, nEntriesGlobal); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(valuesContiguous_); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    LOG(DEBUG) << "\"" << this->name_ << "\" (structured) create valuesContiguous_, nComponents = " << nComponents
      << ", nEntriesLocal = " << nEntriesLocal << ", nEntriesGlobal = " << nEntriesGlobal;
  }

  double *valuesDataContiguous;
  ierr = VecGetArray(valuesContiguous_, &valuesDataContiguous); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    const double *valuesDataComponent;
    ierr = VecGetArrayRead(vectorLocal_[componentNo], &valuesDataComponent); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    VLOG(1) << "  copy " << this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double) << " bytes to contiguous array";
    memcpy(
      valuesDataContiguous + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts(),
      valuesDataComponent,
      this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double)
    );

    ierr = VecRestoreArrayRead(vectorLocal_[componentNo], &valuesDataComponent); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
  }

  ierr = VecRestoreArray(valuesContiguous_, &valuesDataContiguous); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
  valuesContiguousInUse_ = true;
  return valuesContiguous_;
}

//! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
//! this has to be called
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
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
    ierr = VecGetArray(vectorLocal_[componentNo], &valuesDataComponent); CHKERRV(ierr);

    VLOG(1) << "  \"" << this->name_ << "\", component " << componentNo << ", copy " << this->meshPartition_->nDofsLocalWithoutGhosts() << " values, " << this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double) << " bytes from contiguous array";
    memcpy(
      valuesDataComponent,
      valuesDataContiguous + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts(),
      this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double)
    );

    ierr = VecRestoreArray(vectorLocal_[componentNo], &valuesDataComponent); CHKERRV(ierr);
  }

  ierr = VecRestoreArrayRead(valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);
  valuesContiguousInUse_ = false;
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
extractComponent(int componentNo, std::shared_ptr<PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,1>> extractedFieldVariable)
{
  PetscErrorCode ierr;

  // prepare source vector
  Vec vectorSource;
  dof_no_t dofStart = 0;
  // if the contiguous representation is already being used, return contiguous vector
  if (valuesContiguousInUse_)
  {
    vectorSource = valuesContiguous_;
    dofStart = componentNo * this->meshPartition_->nDofsLocalWithoutGhosts();
  }
  else
  {
    vectorSource = vectorLocal_[componentNo];
  }

  double *valuesSource;
  ierr = VecGetArray(vectorSource, &valuesSource); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

  // prepare target vector
  assert(!extractedFieldVariable->valuesContiguousInUse_);

  const double *valuesTarget;
  ierr = VecGetArrayRead(extractedFieldVariable->vectorLocal_[0], &valuesTarget); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

  VLOG(1) << "  copy " << this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double) << " bytes (\"" << this->name_ << "\" component " << componentNo
    << ") to \"" << extractedFieldVariable->name() << "\"";
  memcpy(
    valuesTarget,
    vectorSource + dofStart,
    this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double)
  );

  // restore memory
  ierr = VecRestoreArrayRead(extractedFieldVariable->vectorLocal_[0], &valuesTarget); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
  ierr = VecRestoreArray(vectorSource, &valuesSource); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
}

//! get a vector of local dof nos (from meshPartition), without ghost dofs
template<typename MeshType,typename BasisFunctionType,int nComponents>
std::vector<PetscInt> &PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
localDofNosWithoutGhosts()
{
  assert(this->meshPartition_);
  return this->meshPartition_->localDofNosWithoutGhosts();
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
output(std::ostream &stream)
{
#ifndef NDEBUG  
  // this method gets all local non-ghost values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo, nRanks;
  MPI_Comm_rank(this->meshPartition_->mpiCommunicator(), &ownRankNo);
  MPI_Comm_size(this->meshPartition_->mpiCommunicator(), &nRanks);
  
  // loop over components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    Vec vector = vectorLocal_[componentNo];
    if (valuesContiguous_)
      vector = valuesContiguous_;

    // get global size of vector
    int nEntries, nEntriesLocal;
    PetscErrorCode ierr;
    ierr = VecGetSize(vector, &nEntries); CHKERRV(ierr);
    ierr = VecGetLocalSize(vector, &nEntriesLocal); CHKERRV(ierr);
    
    // retrieve local values
    int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
    std::vector<PetscInt> indices(nDofsLocal);
    if (valuesContiguous_)
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
    VecGetValues(vector, nDofsLocal, indices.data(), localValues.data());
    //VLOG(1) << "localValues: " << localValues;
    
    std::vector<int> localSizes(nRanks);
    localSizes[ownRankNo] = nDofsLocal;
    // MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
    // Note that the recvcount argument at the root indicates the number of items it receives from each process, not the total number of items it receives.
    MPI_Gather(localSizes.data() + ownRankNo, 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, this->meshPartition_->mpiCommunicator());
    
    int maxLocalSize;
    MPI_Allreduce(localSizes.data() + ownRankNo, &maxLocalSize, 1, MPI_INT, MPI_MAX, this->meshPartition_->mpiCommunicator());
    
    //VLOG(1) << "localSizes: " << localSizes << ", maxLocalSize: " << maxLocalSize << ", nRanks: " << nRanks;

    std::vector<double> recvBuffer(maxLocalSize*nRanks);
    std::vector<double> sendBuffer(maxLocalSize,0.0);
    std::copy(localValues.begin(), localValues.end(), sendBuffer.begin());
    
    //VLOG(1) << " sendBuffer: " << sendBuffer;

    MPI_Gather(sendBuffer.data(), maxLocalSize, MPI_DOUBLE, recvBuffer.data(), maxLocalSize, MPI_DOUBLE, 0, this->meshPartition_->mpiCommunicator());
    
    if (ownRankNo == 0)
    {
      if (componentNo == 0)
      {
        stream << "vector \"" << this->name_ << "\" (" << nEntries << " local entries (per component)"
          << (valuesContiguous_? ", valuesContiguous": "") << ")" << std::endl;
      }

      stream << "\"" << this->name_ << "\" component " << componentNo << ": local ordering: [";
      std::vector<double> globalValues(this->meshPartition_->nDofsGlobal());
      for (int rankNo = 0; rankNo < nRanks; rankNo++)
      {
        if (rankNo != 0)
          stream << ",";
        for (dof_no_t dofNoLocal = 0; dofNoLocal < localSizes[rankNo]; dofNoLocal++)
        {
          double value = recvBuffer[rankNo*maxLocalSize + dofNoLocal];
          stream << "  " << value;

        }
      }
      stream << "]," << std::endl;
    }
  }

  if (nRanks > 1)
  {
    // loop over components
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      Vec vector = vectorLocal_[componentNo];
      if (valuesContiguous_)
        vector = valuesContiguous_;

      stream << "locally stored values: [";
      // retrieve local values
      int nDofsLocalWithGhosts = this->meshPartition_->nDofsLocalWithGhosts();
      std::vector<double> localValuesWithGhosts(nDofsLocalWithGhosts);
      std::vector<PetscInt> indices(nDofsLocalWithGhosts);
      for (int i = 0; i < nDofsLocalWithGhosts; i++)
      {
        indices[i] = this->meshPartition_->dofNosLocal()[i] + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts();
      }

      VecGetValues(vector, nDofsLocalWithGhosts, indices.data(), localValuesWithGhosts.data());
      //VLOG(1) << "localValues: " << localValues;

      const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
      for (dof_no_t dofNoLocal = 0; dofNoLocal < this->meshPartition_->nDofsLocalWithoutGhosts(); dofNoLocal++)
      {
        double value = localValuesWithGhosts[dofNoLocal];

        // store value for global dof no
        node_no_t nodeNoLocal = dofNoLocal / nDofsPerNode;
        int dofOnNodeIndex = dofNoLocal % nDofsPerNode;

        std::array<global_no_t,MeshType::dim()> globalCoordinates = this->meshPartition_->getCoordinatesGlobal(nodeNoLocal);
        global_no_t nodeNoGlobal = this->meshPartition_->getNodeNoGlobalNatural(globalCoordinates);

        global_no_t dofNoGlobal = nodeNoGlobal*nDofsPerNode + dofOnNodeIndex;

        stream << "dofNoGlobal=" << dofNoGlobal << ": " << value << ", ";
      }
      stream << "], ghosts: [";

      for (dof_no_t dofNoLocal = this->meshPartition_->nDofsLocalWithoutGhosts(); dofNoLocal < this->meshPartition_->nDofsLocalWithGhosts(); dofNoLocal++)
      {
        double value = localValuesWithGhosts[dofNoLocal];

        // store value for global dof no
        node_no_t nodeNoLocal = dofNoLocal / nDofsPerNode;
        int dofOnNodeIndex = dofNoLocal % nDofsPerNode;

        std::array<global_no_t,MeshType::dim()> globalCoordinates = this->meshPartition_->getCoordinatesGlobal(nodeNoLocal);
        global_no_t nodeNoGlobal = this->meshPartition_->getNodeNoGlobalNatural(globalCoordinates);

        global_no_t dofNoGlobal = nodeNoGlobal*nDofsPerNode + dofOnNodeIndex;

        stream << "dofNoGlobal=" << dofNoGlobal << ": " << value << ", ";
      }
      stream << "]" << std::endl;

      // also output vector using Petsc viewer (not so nice and fails after 1000 files)
#if 0
      PetscViewer viewer;
      static int counter = 0;
      std::stringstream vectorOutputFilename;
      vectorOutputFilename << "vector_" << counter++ << ".txt";
      ierr = PetscViewerASCIIOpen(this->meshPartition_->mpiCommunicator(), vectorOutputFilename.str().c_str(), &viewer); CHKERRV(ierr);
      ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INDEX); CHKERRV(ierr);
      ierr = VecView(vectorGlobal_[componentNo], viewer); CHKERRV(ierr);

      if (ownRankNo == 0)
      {
        stream << "(Vector also written to \"" << vectorOutputFilename.str() << "\".)";
      }
#endif
    }  // componentNo
  }
#endif
}
