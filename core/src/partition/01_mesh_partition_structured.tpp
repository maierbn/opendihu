#include "partition/01_mesh_partition.h"

namespace Partition
{
 
template<typename MeshType,typename BasisFunctionType>
MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
MeshPartition(std::array<global_no_t,MeshType::dim()> globalSize, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), globalSize_(globalSize), hasFullNumberOfNodes_({false})
{
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
 
  PetscErrorCode ierr;
  const int nDofPer1DElement = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  const int ghostLayerWidth = nDofPer1DElement;
  
  // create PETSc DMDA object that is a topology interface handling parallel data layout on structured grids
  if (MeshType::dim() == 1)
  {
    // create 1d decomposition
    ierr = DMDACreate1d(mpiCommunicator(), DM_BOUNDARY_NONE, globalSize_[0], nDofPer1DElement, ghostLayerWidth, 
                        NULL, &dm_); CHKERRV(ierr);
    
    // get global coordinates of local partition
    PetscInt x, m;
    ierr = DMDAGetGhostCorners(dm_, &x, NULL, NULL, &m, NULL, NULL); CHKERRV(ierr);
    beginGlobal_[0] = (global_no_t)x;
    localSizeWithGhosts_[0] = (element_no_t)m;
    
    // get number of ranks in each coordinate direction
    nRanks_[0] = this->rankSubset_->size();
    ierr = DMDAGetInfo(dm_, NULL, NULL, NULL, NULL, &nRanks_[0], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
    
    // get local sizes on the ranks
    const PetscInt *lxData;
    ierr = DMDAGetOwnershipRanges(dm_, &lxData, NULL, NULL);
    
    localSizesOnRanks_[0].assign(lxData,lxData+nRanks_[0]);
    
    // determine if the local partition is at the right border of the global domain
    if (beginGlobal_[0] + localSizeWithGhosts_[0] >= localSize(0))
    {
      hasFullNumberOfNodes_[0] = true;
    }
  }
  else if (MeshType::dim() == 2)
  {
    // create 2d decomposition
    ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                        globalSize_[0], globalSize_[1], PETSC_DECIDE, PETSC_DECIDE,
                        nDofPer1DElement, ghostLayerWidth, NULL, NULL, &dm_); CHKERRV(ierr);
                        
    // get global coordinates of local partition
    PetscInt x, y, m, n;
    ierr = DMDAGetGhostCorners(dm_, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
    beginGlobal_[0] = (global_no_t)x;
    beginGlobal_[1] = (global_no_t)y;
    localSizeWithGhosts_[0] = (element_no_t)m;
    localSizeWithGhosts_[1] = (element_no_t)n;
    
    // get number of ranks in each coordinate direction
    ierr = DMDAGetInfo(dm_, NULL, NULL, NULL, NULL, &nRanks_[0], &nRanks_[1], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
    
    // get local sizes on the ranks
    const PetscInt *lxData;
    const PetscInt *lyData;
    ierr = DMDAGetOwnershipRanges(dm_, &lxData, &lyData, NULL);
    localSizesOnRanks_[0].assign(lxData, lxData + nRanks_[0]);
    localSizesOnRanks_[1].assign(lyData, lyData + nRanks_[1]);
    
    // determine if the local partition is at the x+/y+ border of the global domain
    for (int i = 0; i < MeshType::dim(); i++)
    {
      if (beginGlobal_[i] + localSizeWithGhosts_[i] >= localSize(i))
      {
        hasFullNumberOfNodes_[i] = true;      
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    // create 3d decomposition
    ierr = DMDACreate3d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                        globalSize_[0], globalSize_[1], globalSize_[2],
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                        nDofPer1DElement, ghostLayerWidth, NULL, NULL, NULL, &dm_); CHKERRV(ierr);
                        
    // get global coordinates of local partition
    PetscInt x, y, z, m, n, p;
    ierr = DMDAGetGhostCorners(dm_, &x, &y, &z, &m, &n, &p); CHKERRV(ierr);
    beginGlobal_[0] = (global_no_t)x;
    beginGlobal_[1] = (global_no_t)y;
    beginGlobal_[2] = (global_no_t)z;
    localSizeWithGhosts_[0] = (element_no_t)m;
    localSizeWithGhosts_[1] = (element_no_t)n;
    localSizeWithGhosts_[2] = (element_no_t)p;
    
    // get number of ranks in each coordinate direction
    ierr = DMDAGetInfo(dm_, NULL, NULL, NULL, NULL, &nRanks_[0], &nRanks_[1], &nRanks_[2], NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
    
    // get local sizes on the ranks
    const PetscInt *lxData;
    const PetscInt *lyData;
    const PetscInt *lzData;
    ierr = DMDAGetOwnershipRanges(dm_, &lxData, &lyData, &lzData);
    localSizesOnRanks_[0].assign(lxData, lxData + nRanks_[0]);
    localSizesOnRanks_[1].assign(lyData, lyData + nRanks_[1]);
    localSizesOnRanks_[2].assign(lzData, lzData + nRanks_[2]);
    
    // determine if the local partition is at the x+/y+/z+ border of the global domain
    for (int i = 0; i < MeshType::dim(); i++)
    {
      if (beginGlobal_[i] + localSizeWithGhosts_[i] >= localSize(i))
      {
        hasFullNumberOfNodes_[i] = true;      
      }
    }
  }
  
  this->initializeLocalDofs();
  
  // debugging output
  PetscInt nElements;
  PetscInt nNodes;
  const PetscInt *nodeIndices;
  ierr = DMDAGetElements(dm_, &nElements, &nNodes, &nodeIndices); CHKERRV(ierr);
  
  std::vector<int> nodeIndicesVector(nNodes);
  nodeIndicesVector.assign(nodeIndices, nodeIndices+nNodes);
  
  LOG(DEBUG) << "localSize: " << localSize() << ", globalSize: " << globalSize 
    << ", nElements: " << nElements << ", nNodes: " << nNodes << ", nodeIndices: " << nodeIndicesVector 
    << ", hasFullNumberOfNodes_: " << hasFullNumberOfNodes_;
  LOG(DEBUG) << "nRanks: " << nRanks_ << ", localSizesOnRanks_: " << localSizesOnRanks_ << ", beginGlobal_: " << beginGlobal_ 
    << ", localSizeWithGhosts_: " << localSizeWithGhosts_;
}

template<typename MeshType,typename BasisFunctionType>
MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
MeshPartition(std::array<node_no_t,MeshType::dim()> localSize, std::array<global_no_t,MeshType::dim()> globalSize, 
              std::array<int,MeshType::dim()> beginGlobal, 
              std::array<int,MeshType::dim()> nRanks, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), localSizeWithGhosts_(localSize), globalSize_(globalSize), beginGlobal_(beginGlobal), nRanks_(nRanks)
{
  // partitioning is already prescribed as every rank knows its own local size
 
  // add ghost layer to localSize
  for (int i = 0; i < MeshType::dim(); i++)
  {
    localSizeWithGhosts_[i] += 2;
    beginGlobal_[i] -= 1;
    
    if (beginGlobal_[i] == -1)
    {
      beginGlobal_[i] = 0;
      localSizeWithGhosts_[i]--;
    }
    
    if (beginGlobal_[i] + localSizeWithGhosts_[i] == globalSize_[i])
    {
      localSizeWithGhosts_[i]--;
    }
  }
    
  // determine if the local partition is at the x+/y+ border of the global domain
  for (int i = 0; i < MeshType::dim(); i++)
  {
    if (beginGlobal_[i] + localSizeWithGhosts_[i] >= this->localSize(i))
    {
      hasFullNumberOfNodes_[i] = true;      
    }
  }
 
  // determine localSizesOnRanks_
  std::array<int,MeshType::dim()> ownLocalSize;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    ownLocalSize[i] = this->localSize(i);
    localSizesOnRanks_[i].resize(rankSubset->size());
  }
  
  for (int i = 0; i < MeshType::dim(); i++)
  {
    MPIUtility::handleReturnValue(MPI_Allgather(&ownLocalSize[i], 1, MPI_INT, 
      localSizesOnRanks_[i].data(), rankSubset->size(), MPI_INT, rankSubset->mpiCommunicator()));
  }
  
  LOG(DEBUG) << "localSize: " << this->localSize() << ", globalSize: " << globalSize 
    << ", hasFullNumberOfNodes_: " << hasFullNumberOfNodes_;
  LOG(DEBUG) << "nRanks: " << nRanks_ << ", localSizesOnRanks_: " << localSizesOnRanks_ << ", beginGlobal_: " << beginGlobal_ 
    << ", localSizeWithGhosts_: " << localSizeWithGhosts_;
  
  for (int i = 0; i < MeshType::dim(); i++)
  {
    LOG(DEBUG) << "  beginNodeGlobal(" << i << "): " << beginNodeGlobal(0);
  }
}

//! get the local to global mapping for the current partition
template<typename MeshType,typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
localToGlobalMapping()
{
  ISLocalToGlobalMapping localToGlobalMapping;
  DMGetLocalToGlobalMapping(dm_, &localToGlobalMapping);
  return localToGlobalMapping;
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nRanks(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return nRanks_[coordinateDirection];
}

//! number of entries in the current partition
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
localSize()
{
  element_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    result *= localSizeWithGhosts_[i];
  }
  return result;
}

//! number of entries in the given coordinate direction in the current partition
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
localSize(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return localSizeWithGhosts_[coordinateDirection];
}

//! first local number in current partition
template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
beginGlobal(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return beginGlobal_[coordinateDirection];
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
beginNodeGlobal(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  
  const int nDofPer1DElement = BasisOnMesh::BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
    
  return beginGlobal_[coordinateDirection]*nDofPer1DElement;
}

//! one after last number in current partition
template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
endGlobal(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return beginGlobal_[coordinateDirection] + localSizeWithGhosts_[coordinateDirection];
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
globalSize()
{
  element_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    result *= globalSize_[i];
  }
  return result;
}
  
template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
globalSize(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return globalSize_[coordinateDirection];
}

//! get a vector with the local sizes on every rank, this is needed to create the DMDA in PartitionedPetscVec
template<typename MeshType,typename BasisFunctionType>
std::vector<element_no_t> &MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
localSizesOnRanks(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return localSizesOnRanks_;
}
  
//! get if there are nodes on both borders in the given coordinate direction
//! this is the case if the local partition touches the right/top/back border
template<typename MeshType,typename BasisFunctionType>
bool MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
hasFullNumberOfNodes(int coordinateDirection)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  
  return hasFullNumberOfNodes_[coordinateDirection];
}
  
template<typename MeshType,typename BasisFunctionType>
template<typename T>
void MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
extractLocalNodes(std::vector<T> &vector)
{
  std::vector<T> result(localSize());
  global_no_t resultIndex = 0;
  
  if (MeshType::dim() == 1)
  {
    assert(vector.size() >= beginGlobal_[0] + localSizeWithGhosts_[0]);
    for (global_no_t i = beginGlobal_[0]; i < beginGlobal_[0] + localSizeWithGhosts_[0]; i++)
    {
      result[resultIndex++] = vector[i];
    }
  }
  else if (MeshType::dim() == 2)
  {
    for (global_no_t j = beginGlobal_[1]; j < beginGlobal_[1] + localSizeWithGhosts_[1]; j++)
    {
      for (global_no_t i = beginGlobal_[0]; i < beginGlobal_[0] + localSizeWithGhosts_[0]; i++)
      {
        result[resultIndex++] = vector[j*globalSize(0) + i];
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    for (global_no_t k = beginGlobal_[2]; k < beginGlobal_[2] + localSizeWithGhosts_[2]; k++)
    {
      for (global_no_t j = beginGlobal_[1]; j < beginGlobal_[1] + localSizeWithGhosts_[1]; j++)
      {
        for (global_no_t i = beginGlobal_[0]; i < beginGlobal_[0] + localSizeWithGhosts_[0]; i++)
        {
          result[resultIndex++] = vector[k*globalSize(1)*globalSize(0) + j*globalSize(0) + i];
        }
      }
    }
  }
  
  // store values
  vector.assign(result.begin(), result.end());
}
  
template<typename MeshType,typename BasisFunctionType>
void MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
extractLocalDofs(std::vector<double> &vector)
{
  std::vector<double> result(localSize());
  global_no_t resultIndex = 0;
  
  if (MeshType::dim() == 1)
  {
    assert(vector.size() >= beginGlobal_[0] + localSizeWithGhosts_[0]);
    for (global_no_t i = beginGlobal_[0]; i < beginGlobal_[0] + localSizeWithGhosts_[0]; i++)
    {
      result[resultIndex++] = vector[i];
    }
  }
  else if (MeshType::dim() == 2)
  {
    for (global_no_t j = beginGlobal_[1]; j < beginGlobal_[1] + localSizeWithGhosts_[1]; j++)
    {
      for (global_no_t i = beginGlobal_[0]; i < beginGlobal_[0] + localSizeWithGhosts_[0]; i++)
      {
        result[resultIndex++] = vector[j*globalSize(0) + i];
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    for (global_no_t k = beginGlobal_[2]; k < beginGlobal_[2] + localSizeWithGhosts_[2]; k++)
    {
      for (global_no_t j = beginGlobal_[1]; j < beginGlobal_[1] + localSizeWithGhosts_[1]; j++)
      {
        for (global_no_t i = beginGlobal_[0]; i < beginGlobal_[0] + localSizeWithGhosts_[0]; i++)
        {
          result[resultIndex++] = vector[k*globalSize(1)*globalSize(0) + j*globalSize(0) + i];
        }
      }
    }
  }
  
  // store values
  vector.assign(result.begin(), result.end());
}
  
}  // namespace