#include "partition/partitioned_petsc_vec.h"

//! constructor
template<typename MeshType,typename BasisFunctionType,int nComponents>
PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,MeshType>> meshPartition,
                    std::string name) :
  PartitionedPetscVecBase<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>(meshPartition, name)
{
  //dm_ = meshPartition->dmElements();
  
  createVector();
}

//! constructor, copy from existing vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents> &rhs, std::string name) :
  PartitionedPetscVecBase<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>(rhs.meshPartition_, name)
{
  //dm_ = rhs.dm_;
  
  PetscErrorCode ierr;
  for (int i = 0; i < nComponents; i++)
  {
    ierr = VecCopy(rhs.vectorLocal_[i], vectorLocal_[i]); CHKERRV(ierr);
    ierr = VecCopy(rhs.vectorGlobal_[i], vectorGlobal_[i]); CHKERRV(ierr);
  }
}
  
//! create a distributed Petsc vector, according to partition
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
createVector()
{
  PetscErrorCode ierr;
  
  // The local vector contains the nodal/dof values for the local portion of the current rank. This includes ghost nodes.
  // The global vector manages the whole data and also only stores the local portion of it on the current rank, but without ghost nodes. 
  // If we want to manipulate data, we have to ensure consistency in the parallel global vector (VecGhostUpdateBegin,VecGhostUpdateEnd) 
  // and then fetch the local portion of the global vector together with ghost values into a local vector (VecGhostGetLocalForm),
  // then manipulate the values (using standard VecSetValues/VecGetValues routines) and commit the results (VecGhostRestoreLocalForm, VecGhostUpdateBegin, VecGhostUpdateEnd).
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    const dof_no_t nGhosts = this->meshPartition_->nNodesLocalWithGhosts() - this->meshPartition_->nNodesLocalWithoutGhosts();
    ierr = VecCreateGhost(this->meshPartition_->mpiCommunicator(), this->meshPartition_->nNodesLocalWithoutGhosts(), 
                          this->meshPartition_->nNodesGlobal(), nGhosts, this->meshPartition_->ghostDofGlobalNos().data(), &vectorGlobal_[componentNo]); CHKERRV(ierr);
    
    
    /*ierr = DMCreateGlobalVector(*dm_, &vectorGlobal_[componentNo]); CHKERRV(ierr);
    ierr = DMCreateLocalVector(*dm_, &vectorLocal_[componentNo]); CHKERRV(ierr);*/
    
    // initialize PETSc vector object
    //ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject) vectorGlobal_[componentNo], this->name_.c_str()); CHKERRV(ierr);

    // initialize size of vector
    //ierr = VecSetSizes(vectorLocal_[componentNo], this->meshPartition_->nNodesLocalWithGhosts(), this->meshPartition_->nNodesGlobal()); CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(vectorGlobal_[componentNo]); CHKERRV(ierr);
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
startVectorManipulation()
{
  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;
  
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
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo], &vectorLocal_[componentNo]); CHKERRV(ierr);
    //ierr = DMGlobalToLocalEnd(*dm_, vectorGlobal_[componentNo], INSERT_VALUES, vectorLocal_[componentNo]); CHKERRV(ierr);
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
finishVectorManipulation()
{
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
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecGetValues(vectorLocal_[componentNo], ni, ix, y); CHKERRV(ierr);
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValuesGlobalIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the global vector
  PetscErrorCode ierr;
  ierr = VecGetValues(vectorGlobal_[componentNo], ni, ix, y); CHKERRV(ierr);
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  // this wraps the standard PETSc VecSetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(vectorLocal_[componentNo], ni, ix, y, iora); CHKERRV(ierr);
  
  // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
  // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValue(vectorLocal_[componentNo], row, value, mode); CHKERRV(ierr);
}
/*
//! for a single component vector set all values
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValuesWithGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->meshPartition_->nDofsLocalWithGhosts());
 
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(vectorLocal_[componentNo], values.size(), this->meshPartition_->dofNosLocal().data(), values.data(), petscInsertMode); CHKERRV(ierr);
}

//! for a single component vector set all values, input does not contain values for ghosts
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
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
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
zeroEntries()
{
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecZeroEntries(vectorLocal_[componentNo]); CHKERRV(ierr);
  }
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
Vec &PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
valuesLocal(int componentNo)
{
  return vectorLocal_[componentNo];
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
Vec &PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
valuesGlobal(int componentNo)
{
  return vectorGlobal_[componentNo];
}

//! get a vector of local dof nos (from meshPartition), without ghost dofs
template<typename MeshType,typename BasisFunctionType,int nComponents>
std::vector<PetscInt> &PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
localDofNosWithoutGhosts()
{
  return this->meshPartition_->localDofNosWithoutGhosts();
}