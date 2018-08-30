#include "partition/partitioned_petsc_vec.h"

//! constructor
template<typename MeshType,typename BasisFunctionType,int nComponents>
PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition,
                    std::string name) :
  PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(meshPartition, name), vectorManipulationStarted_(false)
{
  //dm_ = meshPartition->dmElements();
  
  createVector();
}

//! constructor, copy from existing vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(rhs.meshPartition(), name), vectorManipulationStarted_(false)
{
  createVector();
  
  // copy existing values
  PetscErrorCode ierr;
  for (int i = 0; i < nComponents; i++)
  {
    //ierr = VecCopy(rhs.valuesLocal(i), vectorLocal_[i]); CHKERRV(ierr);
    ierr = VecCopy(rhs.valuesGlobal(i), vectorGlobal_[i]); CHKERRV(ierr);
  }
  
  startVectorManipulation();
}
  
//! constructor, copy from existing vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
template<int nComponents2>
PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
PartitionedPetscVec(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents2> &rhs, std::string name) :
  PartitionedPetscVecBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>(rhs.meshPartition(), name), vectorManipulationStarted_(false)
{
  //dm_ = rhs.dm_;
  
  createVector();
  
  // copy existing values
  PetscErrorCode ierr;
  for (int i = 0; i < std::min(nComponents,nComponents2); i++)
  {
    ierr = VecCopy(rhs.valuesGlobal(i), vectorGlobal_[i]); CHKERRV(ierr);
  }
  
  vectorManipulationStarted_ = false;
  startVectorManipulation();
}
  
//! create a distributed Petsc vector, according to partition
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
createVector()
{
  VLOG(2) << "\"" << this->name_ << "\" createVector with " << nComponents << " components, size local: " << this->meshPartition_->nNodesLocalWithoutGhosts() 
    << ", global: " << this->meshPartition_->nNodesGlobal() << ", ghost dof global nos: " << this->meshPartition_->ghostDofGlobalNos();
  PetscErrorCode ierr;
  
  // The local vector contains the nodal/dof values for the local portion of the current rank. This includes ghost nodes.
  // The global vector manages the whole data and also only stores the local portion of it on the current rank, but without ghost nodes. 
  // If we want to manipulate data, we have to ensure consistency in the parallel global vector (VecGhostUpdateBegin,VecGhostUpdateEnd) 
  // and then fetch the local portion of the global vector together with ghost values into a local vector (VecGhostGetLocalForm),
  // then manipulate the values (using standard VecSetValues/VecGetValues routines) and commit the results (VecGhostRestoreLocalForm, VecGhostUpdateBegin, VecGhostUpdateEnd).
  
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    const dof_no_t nGhostDofs = this->meshPartition_->nDofsLocalWithGhosts() - this->meshPartition_->nDofsLocalWithoutGhosts();
    ierr = VecCreateGhost(this->meshPartition_->mpiCommunicator(), this->meshPartition_->nDofsLocalWithoutGhosts(),
                          this->meshPartition_->nDofsGlobal(), nGhostDofs, this->meshPartition_->ghostDofGlobalNos().data(), &vectorGlobal_[componentNo]); CHKERRV(ierr);
    
    
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
  
  startVectorManipulation();
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
startVectorManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" startVectorManipulation";
  
  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;
  if (vectorManipulationStarted_)
  {
    LOG(DEBUG) << "\"" << this->name_ << "\", startVectorManipulation called multiple times without finishVectorManipulation in between. (Creation of the vector also calls startVectorManipulation).";
  }
  
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
  vectorManipulationStarted_ = true;
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
finishVectorManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" finishVectorManipulation";
  
  if (!vectorManipulationStarted_)
  {
    LOG(ERROR) << "\"" << this->name_ << "\", finishVectorManipulation called without previous startVectorManipulation";
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
  vectorManipulationStarted_ = false;
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecGetValues(vectorLocal_[componentNo], ni, ix, y); CHKERRV(ierr);
  
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

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
getValuesGlobalIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the global vector
  PetscErrorCode ierr;
  ierr = VecGetValues(vectorGlobal_[componentNo], ni, ix, y); CHKERRV(ierr);
  
  // debugging output
  if (VLOG_IS_ON(3))
  {
    std::stringstream str;
    str << "\"" << this->name_ << "\" getValuesGlobalIndexing(componentNo=" << componentNo << ", indices=";
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
    str << (iora == INSERT_VALUES? "INSERT_VALUES" : "ADD_VALUES");
    str << ")";
    VLOG(3) << str.str();
  }
  
  // this wraps the standard PETSc VecSetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(vectorLocal_[componentNo], ni, ix, y, iora); CHKERRV(ierr);
  
  // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
  // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
}

template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  VLOG(3) << "\"" << this->name_ << "\" setValue(componentNo=" << componentNo << ", row=" << row << ", value=" << value
    << (mode == INSERT_VALUES? "INSERT_VALUES" : "ADD_VALUES") << ")";
  
  // this wraps the standard PETSc VecSetValue on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValue(vectorLocal_[componentNo], row, value, mode); CHKERRV(ierr);
}

//! set values from another vector
template<typename MeshType,typename BasisFunctionType,int nComponents>
void PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents,Mesh::isStructured<MeshType>>::
setValues(PartitionedPetscVec<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,nComponents> &rhs)
{
  VLOG(3) << "\"" << this->name_ << "\" setValues(rhs \"" << rhs.name() << "\"), this calls startVectorManipulation()";
  
  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    ierr = VecCopy(rhs.valuesGlobal(componentNo), vectorGlobal_[componentNo]);
  }
  
  startVectorManipulation();
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
  return vectorGlobal_[componentNo];
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
  // this method gets all values and outputs them to stream, only on rank 0
  PetscMPIInt ownRankNo, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
  MPI_Comm_size(PETSC_COMM_WORLD, &nRanks);
  
  startVectorManipulation();

  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // get global size of vector
    int nEntries, nEntriesLocal;
    PetscErrorCode ierr;
    ierr = VecGetSize(vectorLocal_[componentNo], &nEntries); CHKERRV(ierr);
    ierr = VecGetLocalSize(vectorLocal_[componentNo], &nEntriesLocal); CHKERRV(ierr);
    
    // retrieve local values
    int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
    std::vector<double> localValues(nDofsLocal);
    
    VecGetValues(vectorLocal_[componentNo], nDofsLocal, this->meshPartition_->dofNosLocal().data(), localValues.data());
    //VLOG(1) << "localValues: " << localValues;
    
    std::vector<int> localSizes(nRanks);
    localSizes[ownRankNo] = nDofsLocal;
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
        stream << "vector \"" << this->name_ << "\" (" << nEntries << " global entries per component)" << std::endl;
      }

      stream << "\"" << this->name_ << "\" component " << componentNo << ": [";
      for (int rankNo = 0; rankNo < nRanks; rankNo++)
      {
        if (rankNo != 0)
          stream << ",";
        for (int i = 0; i < localSizes[rankNo]; i++)
        {
          stream << "  " << recvBuffer[rankNo*maxLocalSize + i];
        }
      }
      stream << "]" << std::endl;
    }
  }  // componentNo

  finishVectorManipulation();
#endif
}
