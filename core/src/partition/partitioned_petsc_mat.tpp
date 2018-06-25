#include "partition/partitioned_petsc_vec.h"

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
PartitionedPetscMat(std::shared_ptr<MeshPartition> meshPartition, int nComponents) :
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
}

//! constructor
template<int D, typename MeshType, typename BasisFunctionType>
PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
PartitionedPetscMat(Mat &matrix) :
  meshPartition_(nullptr), nComponents_(-1), matrix_(matrix)
{
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

    ierr = MatCreate(rankSubset_->mpiCommunicator(), &matrix_); CHKERRV(ierr);
    ierr = MatSetSizes(matrix_, partition_.localSize(), partition_.localSize(), meshPartition->globalSize(), meshPartition->globalSize()); CHKERRV(ierr);
    ierr = MatSetFromOptions(matrix_); CHKERRV(ierr);                        
    
    // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
    //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

    // dense matrix
    //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);

    // sparse matrix
    ierr = MatMPIAIJSetPreallocation(matrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
    ierr = MatSeqAIJSetPreallocation(matrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
  }
  else 
  {
    // parallel API
    ierr = DMSetMatrixPreallocateOnly(dm_, true); CHKERRV(ierr);  // do not fill zero entries when DMCreateMatrix is called
    ierr = DMCreateMatrix(dm_, &matrix_); CHKERRV(ierr);
  }
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
matSetValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
{
  // this wraps the standard PETSc MatSetValue on the local matrix
  PetscErrorCode ierr;
  
  ierr = MatSetValuesLocal(matrix_, 1, &row, 1, &col, &value, mode); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
matSetValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
{
  // this wraps the standard PETSc MatSetValues on the local matrix
  PetscErrorCode ierr;
  
  ierr = MatSetValuesLocal(matrix_, m, idxm, n, idxn, v, addv); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
matZeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag)
{
  PetscErrorCode ierr;
  ierr = MatZeroRowsColumnsLocal(matrix_, numRows, rows, diag, NULL, NULL); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
matAssembly(MatAssemblyType type)
{
  // this wraps the standard PETSc assembleBegin/End
  PetscErrorCode ierr;
  
  ierr = MatAssemblyBegin(matrix_, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrix_, type); CHKERRV(ierr);
}

template<int D, typename MeshType, typename BasisFunctionType>
void PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
matGetValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[])
{
  // this wraps the standard PETSc MatGetValues, for the global indexing, only retrieves locally stored indices
  PetscErrorCode ierr;
  
  ierr = MatGetValues(matrix_, m, idxm, n, idxn, v); CHKERRV(ierr);
}


//! get a reference to the PETSc matrix
template<int D, typename MeshType, typename BasisFunctionType>
Mat &PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
values()
{
  return this->matrix_;
}