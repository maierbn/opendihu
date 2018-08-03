#include "partition/01_mesh_partition.h"

#include <petscsys.h>

#include "mesh/mesh.h"
#include "utility/mpi_utility.h"

namespace Partition
{

MeshPartition<Mesh::None>::
MeshPartition(global_no_t globalSize, std::shared_ptr<RankSubset> rankSubset) :
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
    localSize_++;
  
  beginGlobal_ = rankNo * sizePerRank + std::min(residual, (global_no_t)rankNo);
  
  // initialize local dofs list 
  this->initializeLocalDofs();
}

//! get the local to global mapping for the current partition
ISLocalToGlobalMapping MeshPartition<Mesh::None>::
localToGlobalMapping()
{
  PetscErrorCode ierr;
  std::vector<PetscInt> globalDofNos(localSize());
  std::iota(globalDofNos.begin(), globalDofNos.end(), beginGlobal_);
  ISLocalToGlobalMapping localToGlobalMapping;
  ierr = ISLocalToGlobalMappingCreate(mpiCommunicator(), 1, localSize(), 
                                      globalDofNos.data(), PETSC_COPY_VALUES, &localToGlobalMapping); CHKERRABORT(mpiCommunicator(), ierr);

  return localToGlobalMapping;
}
 
//! number of entries in the current partition
element_no_t MeshPartition<Mesh::None>::
localSize()
{
  return localSize_;
}

//! number of entries in total
global_no_t MeshPartition<Mesh::None>::
globalSize()
{
  return globalSize_;
}
    
void MeshPartition<Mesh::None>::
extractLocalDofs(std::vector<double> &vector)
{
  extractLocalNodes<double>(vector);
}
  
}  // namespace