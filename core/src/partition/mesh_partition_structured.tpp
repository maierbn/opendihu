#include "partition/mesh_partition.h"

namespace Partition
{
 
template<int D, typename MeshType, typename BasisFunctionType>
MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
MeshPartition(std::array<node_no_t,D> globalSize, std::shared_ptr<RankSubset> rankSubset) :
  globalSize_(globalSize), rankSubset_(rankSubset)
{
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
 
  PetscErrorCode ierr;
  
  // create PETSc DMDA object that is a topology interface handling parallel data layout on structured grids
  if (D == 1)
  {
    // create 1d decomposition
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate1d(mpiCommunicator(), DM_BOUNDARY_NONE, globalSize_[0], nComponents_, ghostLayerWidth, 
                        NULL, &dm_); CHKERRV(ierr);
    
    // get global coordinates of local partition
    Petsc_Int x, m;
    ierr = DMDAGetGhostCorners(dm_, &x, NULL, NULL, &m, NULL, NULL); CHKERRV(ierr);
    beginGlobal_[0] = (global_no_t)x;
    localSizeWithGhosts_[0] = (element_no_t)m;
    
    // get number of ranks in each coordinate direction
    nRanks_[0] = rankSubset_->size();
    ierr = DMDAGetInfo(dm_, NULL, NULL, NULL, NULL, &nRanks_[0], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
    
    localSizesOnRanks_[0].resize(nRanks_[0]);
    
    // get local sizes on the ranks
    std::vector<PetscInt> lx;
    ierr = DMDAGetOwnershipRanges(dm_, lx.data(), NULL, NULL);
    localSizesOnRanks_[0] = lx;
  }
  else if (D == 2)
  {
    // create 2d decomposition
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                        globalSize_[0], globalSize_[1], PETSC_DECIDE, PETSC_DECIDE,
                        nComponents_, ghostLayerWidth, NULL, NULL, &dm_); CHKERRV(ierr);
                        
    // get global coordinates of local partition
    Petsc_Int x, y, m, n;
    ierr = DMDAGetGhostCorners(dm_, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
    beginGlobal_[0] = (global_no_t)x;
    beginGlobal_[1] = (global_no_t)y;
    localSizeWithGhosts_[0] = (element_no_t)m;
    localSizeWithGhosts_[1] = (element_no_t)n;
    
    // get number of ranks in each coordinate direction
    ierr = DMDAGetInfo(dm_, NULL, NULL, NULL, NULL, &nRanks_[0], &nRanks_[1], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
    
    localSizesOnRanks_[0].resize(nRanks_[0]);
    localSizesOnRanks_[1].resize(nRanks_[1]);
    
    // get local sizes on the ranks
    std::vector<PetscInt> lx;
    std::vector<PetscInt> ly;
    ierr = DMDAGetOwnershipRanges(dm_, lx.data(), ly.data(), NULL);
    localSizesOnRanks_[0] = lx;
    localSizesOnRanks_[1] = ly;
  }
  else if (D == 3)
  {
    // create 3d decomposition
    int ghostLayerWidth = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::nAverageNodesPerElement();
    ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                        globalSize_[0], globalSize_[1], globalSize_[2],
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                        nComponents_, ghostLayerWidth, NULL, NULL, NULL, &dm_); CHKERRV(ierr);
                        
    // get global coordinates of local partition
    Petsc_Int x, y, z, m, n, p;
    ierr = DMDAGetGhostCorners(dm_, &x, &y, &z, &m, &n, &p); CHKERRV(ierr);
    beginGlobal_[0] = (global_no_t)x;
    beginGlobal_[1] = (global_no_t)y;
    beginGlobal_[2] = (global_no_t)z;
    localSizeWithGhosts_[0] = (element_no_t)m;
    localSizeWithGhosts_[1] = (element_no_t)n;
    localSizeWithGhosts_[2] = (element_no_t)p;
    
    // get number of ranks in each coordinate direction
    ierr = DMDAGetInfo(dm_, NULL, NULL, NULL, NULL, &nRanks_[0], &nRanks_[1], &nRanks_[2], NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
    
    localSizesOnRanks_[0].resize(nRanks_[0]);
    localSizesOnRanks_[1].resize(nRanks_[1]);
    localSizesOnRanks_[2].resize(nRanks_[2]);
    
    // get local sizes on the ranks
    std::vector<PetscInt> lx;
    std::vector<PetscInt> ly;
    std::vector<PetscInt> lz;
    ierr = DMDAGetOwnershipRanges(dm_, lx.data(), ly.data(), lz.data());
    localSizesOnRanks_[0] = lx;
    localSizesOnRanks_[1] = ly;
    localSizesOnRanks_[2] = lz;
  }
  
}

template<int D, typename MeshType, typename BasisFunctionType>
AO &MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
applicationOrdering()
{
  PetscErrorCode ierr;
 
  // retrieve the local to global mapping
  ierr = DMDAGetAO(dm_, &ao_); CHKERRV(ierr);
  
  if (VLOG_IS_ON(1))
  {
    PetscViewer viewer;
    ierr = PetscViewerCreate(mpiCommunicator(),PetscViewer &viewer); CHKERRV(ierr);
    VLOG(1) << "Application ordering: ";
    ierr = AOView(ao_); CHKERRV(ierr);
  }
  return this->ao_;
}

template<int D, typename MeshType, typename BasisFunctionType>
int MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
nRanks(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < D);
  return nRanks_[coordinateDirection];
}

//! number of entries in the current partition
template<int D, typename MeshType, typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
localSize()
{
  element_no_t result = 1;
  for (int i = 0; i < D; i++)
  {
    localSize *= localSizeWithGhosts_[i];
  }
  return result;
}

//! number of entries in the given coordinate direction in the current partition
template<int D, typename MeshType, typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
localSize(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < D);
  return localSizeWithGhosts_[coordinateDirection];
}

//! first local number in current partition
template<int D, typename MeshType, typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
beginGlobal(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < D);
  return beginGlobal_[coordinateDirection];
}

//! one after last number in current partition
template<int D, typename MeshType, typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
endGlobal(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < D);
  return beginGlobal_[coordinateDirection] + localSizeWithGhosts_[coordinateDirection];
}
  
  globalSize_
  
template<int D, typename MeshType, typename BasisFunctionType>
MPI_Comm MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
mpiCommunicator()
{
  return rankSubset_->mpiCommunicator();
}

template<int D, typename MeshType, typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
global_no_t globalSize()
{
  element_no_t result = 1;
  for (int i = 0; i < D; i++)
  {
    result *= globalSize_[i];
  }
  return result;
}
  
template<int D, typename MeshType, typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
global_no_t globalSize(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < D);
  return globalSize_[coordinateDirection];
}

//! get a vector with the local sizes on every rank
std::vector<element_no_t> &MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
localSizesOnRanks(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < D);
  return localSizesOnRanks_;
}
  
template <typename T>
void MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>::
extractLocalNumbers(std::vector<T> &vector)
{
  std::vector<T> result(localSize());
  global_no_t resultIndex = 0;
  
  if (D == 1)
  {
    for (global_no_t i = beginGlobal_[0]; i < localSizeWithGhosts_[0]; i++)
    {
      result[resultIndex++] = vector[i];
    }
  }
  else if (D == 2)
  {
    for (global_no_t j = beginGlobal_[1]; j < localSizeWithGhosts_[1]; j++)
    {
      for (global_no_t i = beginGlobal_[0]; i < localSizeWithGhosts_[0]; i++)
      {
        result[resultIndex++] = vector[i];
      }
    }
  }
  else if (D == 3)
  {
    for (global_no_t k = beginGlobal_[2]; k < localSizeWithGhosts_[2]; k++)
    {
      for (global_no_t j = beginGlobal_[1]; j < localSizeWithGhosts_[1]; j++)
      {
        for (global_no_t i = beginGlobal_[0]; i < localSizeWithGhosts_[0]; i++)
        {
          result[resultIndex++] = vector[k*globalSize(2)...];
        }
      }
    }
  }
  
   = (global_no_t)x;
    beginGlobal_[1] = (global_no_t)y;
    beginGlobal_[2] = (global_no_t)z;
    localSizeWithGhosts_[0] = (element_no_t)m;
    localSizeWithGhosts_[1] = (element_no_t)n;
    localSizeWithGhosts_[2] = (element_no_t)p;
    
}
  
}  // namespace