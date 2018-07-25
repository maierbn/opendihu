#include "partition/01_mesh_partition.h"

#include "mesh/mesh.h"

namespace Partition
{
 
template<int D, typename BasisFunctionType>
MeshPartition<Mesh::None>::
MeshPartition(element_no_t globalSize, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), globalSize_(globalSize)
{
  
  
  // get the own current MPI rank
  int rankNo;
  int nRanks;
  MPIUtility::handleReturnValue (MPI_Comm_size(this->mpiCommunicator(), &nRanks));
  MPIUtility::handleReturnValue (MPI_Comm_rank(this->mpiCommunicator(), &rankNo));
  
  global_no_t sizePerRank = globalSize_ / nRanks;
  global_no_t residual = globalSize_ - sizePerRank*nRanks;
  localSize_ = element_no_t(sizePerRank);
  
  if (rankNo < residual)
    localSize++;
  
  beginGlobal_ = rankNo * sizePerRank + std::min(residual, (global_no_t)rankNo);
  
  // initialize local dofs list 
  this->initializeLocalDofs();
}

template<int D, typename BasisFunctionType>
AO &MeshPartition<Mesh::None>::
applicationOrdering()
{
}

//! get the local to global mapping for the current partition
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<Mesh::None>::
localToGlobalMapping()
{
  PetscErrorCode ierr;
  std::vector<PetscInt> globalDofNos(localSize());
  std::iota(globalDofNos.begin(), globalDofNos.end(), beginGlobal_);
  ISLocalToGlobalMapping localToGlobalMapping;
  ierr = ISLocalToGlobalMappingCreate(mpiCommunicator(), 1, localSize(), 
                                      globalDofNos.size(), PETSC_COPY_VALUES, &localToGlobalMapping); CHKERRABORT(ierr);

  return localToGlobalMapping;
}
 
//! number of entries in the current partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<Mesh::None>::
localSize()
{
  return localSize_;
}

//! number of entries in total
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<Mesh::None>::
globalSize()
{
  return globalSize_;
}
  
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<Mesh::None>::
extractLocalNumbers(std::vector<T> &vector)
{
  // copy the entries at [beginGlobal_, beginGlobal_+localSize_] to [0,localSize_]
  for (int i = 0; i < localSize_; i++)
  {
    vector[i] = vector[i+beginGlobal_];
  }
  
  // crop the vector at localSize_
  vector.resize(localSize_);
}
  
}  // namespace