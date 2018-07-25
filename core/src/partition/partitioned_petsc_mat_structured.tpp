#include "partition/partitioned_petsc_vec.h"

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
PartitionedPetscMat(std::shared_ptr<MeshPartition> meshPartition, int nComponents,
                    int diagonalNonZeros, int offdiagonalNonZeros) :
  meshPartition_(meshPartition), nComponents_(nComponents)
{
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
  
  PetscErrorCode ierr;
  
  // create PETSc DMDA object that is a topology interface handling parallel data layout on structured grids
  // This also contains the number of components for each dof. Therefore we can't simply use the DM object of meshPartition, but have to create a new DM object here.
  if (D == 1)
  {
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate1d(meshPartition_->mpiCommunicator(), DM_BOUNDARY_NONE, meshPartition->globalSize(0), nComponents_, ghostLayerWidth, 
                        NULL, dm_); CHKERRV(ierr);
  }
  else if (D == 2)
  {
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate2d(meshPartition_->mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                        meshPartition->globalSize(0), meshPartition->globalSize(1), PETSC_DECIDE, PETSC_DECIDE,
                        nComponents_, ghostLayerWidth, NULL, NULL, dm_); CHKERRV(ierr);
  }
  else if (D == 3)
  {
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate2d(meshPartition_->mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                        meshPartition->globalSize(0), meshPartition->globalSize(1), meshPartition->globalSize(2), 
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                        nComponents_, ghostLayerWidth, NULL, NULL, NULL, dm_); CHKERRV(ierr);
  }
  
  createMatrix(diagonalNonZeros, offdiagonalNonZeros);
}

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
PartitionedPetscMat(Mat &matrix) :
  meshPartition_(nullptr), nComponents_(-1), matrix_(matrix)
{ 
  //! constructor to simply wrap an existing Mat, as needed in nonlinear solver callback functions for jacobians
}

//! create a distributed Petsc matrix, according to the given partition
template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
createMatrix(int diagonalNonZeros, int offdiagonalNonZeros)
{
  PetscErrorCode ierr;
  
  const bool serial = true;   /// use the serial PETSc API
  
  // sparse matrix, serial
  if (serial)
  {
    //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.localSize(), partition.localSize(), n, n,
    //                    diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &matrix); CHKERRV(ierr);

    assert(meshPartition_);
   
    ierr = MatCreate(meshPartition_->mpiCommunicator(), &matrix_); CHKERRV(ierr);
    ierr = MatSetSizes(matrix_, meshPartition_->localSize(), meshPartition_->localSize(), meshPartition->globalSize(), meshPartition->globalSize()); CHKERRV(ierr);
    ierr = MatSetFromOptions(matrix_); CHKERRV(ierr);                        
    
    // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
    //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

    // dense matrix
    //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

    // sparse matrix
    ierr = MatMPIAIJSetPreallocation(matrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
    ierr = MatSeqAIJSetPreallocation(matrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
    
    ierr = MatSetLocalToGlobalMapping(matrix_, meshPartition_->localToGlobalMapping(), meshPartition_->localToGlobalMapping()); CHKERRV(ierr);
  }
  else 
  {
    // parallel API
    ierr = DMSetMatrixPreallocateOnly(dm_, true); CHKERRV(ierr);  // do not fill zero entries when DMCreateMatrix is called
    ierr = DMCreateMatrix(dm_, &matrix_); CHKERRV(ierr);
  }
}