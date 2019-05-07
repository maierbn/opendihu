#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

//! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
//! after manipulation of the vector has finished one has to call restoreValuesContiguous
template<typename MeshType,typename BasisFunctionType,int nComponents>
Vec &PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValuesContiguous()
{
  VLOG(2) << "\"" << this->name_ << "\" getValuesContiguous(), nComponents=" << nComponents;

  this->setRepresentationContiguous();
  return this->valuesContiguous_;
}

//! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
//! after manipulation of the vector has finished one has to call restoreValuesContiguous
template<typename MeshType,typename BasisFunctionType>
Vec &PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,1,Mesh::isStructured<MeshType>>::
getValuesContiguous()
{
  VLOG(2) << "\"" << this->name_ << "\" getValuesContiguous(), nComponents=1";

  return this->vectorGlobal_[0];
}


//! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVecNComponentsStructured where there is one local vector with ghosts for each component.
//! this has to be called
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
restoreValuesContiguous()
{
  VLOG(2) << "\"" << this->name_ << "\" restoreValuesContiguous() nComponents=" << nComponents;

  assert(this->valuesContiguous_ != PETSC_NULL);
  if (this->currentRepresentation_ != Partition::values_representation_t::representationContiguous)
  {
    LOG(FATAL) << "Called restoreValuesContiguous() in representation "
      << Partition::valuesRepresentationString[this->currentRepresentation_] << ", probably without previous getValuesContiguous()";
  }

  PetscErrorCode ierr;
  const double *valuesDataContiguous;
  ierr = VecGetArrayRead(this->valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    double *valuesDataComponent;
    ierr = VecGetArray(this->vectorLocal_[componentNo], &valuesDataComponent); CHKERRV(ierr);

    VLOG(1) << "  \"" << this->name_ << "\", component " << componentNo << ", copy " << this->meshPartition_->nDofsLocalWithoutGhosts() << " values, " << this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double) << " bytes from contiguous array";
    memcpy(
      valuesDataComponent,
      valuesDataContiguous + componentNo*this->meshPartition_->nDofsLocalWithoutGhosts(),
      this->meshPartition_->nDofsLocalWithoutGhosts()*sizeof(double)
    );

    ierr = VecRestoreArray(this->vectorLocal_[componentNo], &valuesDataComponent); CHKERRV(ierr);
  }

  ierr = VecRestoreArrayRead(this->valuesContiguous_, &valuesDataContiguous); CHKERRV(ierr);
  this->currentRepresentation_ = Partition::values_representation_t::representationLocal;
}

//! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVecNComponentsStructured where there is one local vector with ghosts for each component.
//! this has to be called
template<typename MeshType,typename BasisFunctionType>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,1,Mesh::isStructured<MeshType>>::
restoreValuesContiguous()
{
  VLOG(2) << "\"" << this->name_ << "\" restoreValuesContiguous() nComponents=1";
  // if there is only one component, do not use the contiguous vector
}
