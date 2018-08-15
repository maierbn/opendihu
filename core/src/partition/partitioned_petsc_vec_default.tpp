#include "partition/partitioned_petsc_vec.h"

template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
PartitionedPetscVec<BasisOnMeshType,nComponents,DummyForTraits>::
PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition, std::string name) :
  PartitionedPetscVecBase<BasisOnMeshType>(meshPartition)
{
  PetscErrorCode ierr;
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // initialize PETSc vector object
    ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &values_[componentNo]); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject) values_[componentNo], name.c_str()); CHKERRV(ierr);

    // initialize size of vector
    ierr = VecSetSizes(values_[componentNo], this->meshPartition_->localSize(), this->meshPartition_->globalSize()); CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(values_[componentNo]); CHKERRV(ierr);
  }
}

//! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
startVectorManipulation()
{
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
finishVectorManipulation()
{
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

//! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  // this wraps the standard PETSc VecSetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(values_[componentNo], ni, ix, y, iora); CHKERRV(ierr);
}

//! wrapper to the PETSc VecSetValue, acting only on the local data
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValue(values_[componentNo], row, value, mode); CHKERRV(ierr);
}

//! for a single component vector set all values
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
setValues(int componentNo, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->meshPartition_->localSize());
 
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValue(values_[componentNo], this->meshPartition_->localDofs().data(), values.data(), petscInsertMode); CHKERRV(ierr);
}
  
//! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecGetValues(values_[componentNo], ni, ix, y); CHKERRV(ierr); 
}

//! wrapper to the PETSc VecGetValues, on the global vector with global indexing
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
getValuesGlobalIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  getValues(componentNo, ni, ix, y);
}

//! get all locally stored values
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
getLocalValues(int componentNo, std::vector<double> &values)
{
  VecGetValues(values_[componentNo], this->meshPartition_->localSize(), this->meshPartition_->localDofs().data(), values.data());
}

//! set all entries to zero, wraps VecZeroEntries
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
void PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
zeroEntries()
{
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecZeroEntries(values_[componentNo]); CHKERRV(ierr);
  }
}

//! get the local Vector of a specified component
template<typename BasisOnMeshType, int nComponents, typename DummyForTraits>
Vec &PartitionedPetscVec<BasisOnMeshType, nComponents, DummyForTraits>::
values(int componentNo)
{
  return values_[componentNo];
}
