#include "partition/partitioned_petsc_vec.h"

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
PartitionedPetscVec(MeshPartition &meshPartition, int nComponents, std::string name) :
  meshPartition_(meshPartition), nComponents_(nComponents)
{
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
 
  PetscErrorCode ierr;
  
  // Create PETSc DMDA object (distributed array). 
  // This is a topology interface handling parallel data layout on structured grids.
  if (D == 1)
  {
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate1d(meshPartition_->mpiCommunicator(), DM_BOUNDARY_NONE, meshPartition->globalSize(0), nComponents_, ghostLayerWidth, 
                        meshPartition.localSizesOnRanks(0).data(), dm_); CHKERRV(ierr);
  }
  else if (D == 2)
  {
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate2d(meshPartition_->mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                 meshPartition->globalSize(0), meshPartition->globalSize(1), meshPartition->nRanks(0), meshPartition->nRanks(1),
                 nComponents_, ghostLayerWidth, 
                 meshPartition.localSizesOnRanks(0).data(), meshPartition.localSizesOnRanks(1).data(),
                 dm_); CHKERRV(ierr);
  }
  else if (D == 3)
  {
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate2d(meshPartition_->mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                 meshPartition->globalSize(0), meshPartition->globalSize(1), meshPartition->globalSize(2), 
                 meshPartition->nRanks(0), meshPartition->nRanks(1), meshPartition->nRanks(2),
                 nComponents_, ghostLayerWidth, meshPartition.localSizesOnRanks(0).data(), meshPartition.localSizesOnRanks(1).data(),
                 meshPartition.localSizesOnRanks(2).data(), dm_); CHKERRV(ierr);
  }
  
  createVector(name);
}

//! create a distributed Petsc vector, according to partition
template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
createVector(std::string name)
{
  PetscErrorCode ierr;
  
  // The local vector contains the nodal/dof values for the local portion of the current rank. This includes ghost nodes, the correct amount of variables is assured by the DM object.
  // The global vector manages the whole data and also only stores the local portion of it on the current rank, but without ghost nodes. 
  // If we want to manipulate data, we have to fetch the global data into the local vector (DMGlobalToLocalBegin/DMGlobalToLocalEnd),
  // then manipulate the values (using standard VecSetValues/VecGetValues routines) and commit the results (DMLocalToGlobalBegin/DMLocalToGlobalEnd).
  
  ierr = DMCreateGlobalVector(dm_, &vectorGlobal_); CHKERRV(ierr);
  ierr = DMCreateLocalVector(dm_, &vectorLocal_); CHKERRV(ierr);
  
  // serial vector
  if (false)
  {
    // initialize PETSc vector object
    ierr = VecCreate(rankSubset_->mpiCommunicator(), &localVector); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject) localVector, name.c_str()); CHKERRV(ierr);

    // initialize size of vector
    ierr = VecSetSizes(localVector, partition.localSize(), nEntries); CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(localVector); CHKERRV(ierr);
  }
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
startVectorManipulation()
{
  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;
  ierr = VecZeroEntries(vectorLocal_); CHKERRV(ierr);
  ierr = DMGlobalToLocalBegin(dm_, vectorGlobal_, INSERT_VALUES, vectorLocal_); CHKERRV(ierr);
  ierr = DMGlobalToLocalEnd(dm_, vectorGlobal_, INSERT_VALUES, vectorLocal_); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
finishVectorManipulation()
{
  // Copy the local values vectors into the global vector. ADD_VALUES means that ghost values are reduced (summed up)
  PetscErrorCode ierr;
  ierr = VecZeroEntries(vectorGlobal_); CHKERRV(ierr);
  ierr = DMLocalToGlobalBegin(dm_, vectorLocal_, ADD_VALUES, vectorGlobal_); CHKERRV(ierr);
  ierr = DMLocalToGlobalEnd(dm_, vectorLocal_, ADD_VALUES, vectorGlobal_); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
getValues(PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecGetValues(vectorLocal_, ni, ix, y); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
getValuesGlobalIndexing(PetscInt ni, const PetscInt ix[], PetscScalar y[])
{
  // this wraps the standard PETSc VecGetValues on the global vector
  PetscErrorCode ierr;
  ierr = VecGetValues(vectorGlobal_, ni, ix, y); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
setValues(PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  // this wraps the standard PETSc VecSetValues on the local vector
  PetscErrorCode ierr;
  ierr = VecSetValues(vectorLocal_, ni, ix, y, iora); CHKERRV(ierr);
  
  // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
  // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
zeroEntries()
{
  PetscErrorCode ierr;
  ierr = VecZeroEntries(vectorLocal_); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
Vec &PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
values()
{
  return vectorLocal_;
}