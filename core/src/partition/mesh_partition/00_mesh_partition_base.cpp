#include "partition/mesh_partition/00_mesh_partition_base.h"

#include <numeric>

#include "utility/string_utility.h"
#include "utility/vector_operators.h"
#include "easylogging++.h"


namespace Partition
{
 
MeshPartitionBase::MeshPartitionBase(std::shared_ptr<RankSubset> rankSubset) :
  rankSubset_(rankSubset)
{
}

MeshPartitionBase::~MeshPartitionBase()
{
}

int MeshPartitionBase::nRanks() const
{
  return this->rankSubset_->size();
}

int MeshPartitionBase::ownRankNo()
{
  return this->rankSubset_->ownRankNo();
}

std::shared_ptr<RankSubset> MeshPartitionBase::rankSubset() const
{
  return rankSubset_;
}

MPI_Comm MeshPartitionBase::mpiCommunicator() const
{
  return rankSubset_->mpiCommunicator();
}

//! fill the dofNosLocal_ vectors
void MeshPartitionBase::createLocalDofOrderings()
{
  dof_no_t nDofsLocalWithGhosts = this->nDofsLocalWithGhosts();
  dof_no_t nDofsLocalWithoutGhosts = this->nDofsLocalWithoutGhosts();

  // fill dofNosLocal_ 
  dofNosLocal_.resize(nDofsLocalWithGhosts);
  std::iota(dofNosLocal_.begin(), dofNosLocal_.end(), 0);
  
  LOG(DEBUG) << " MeshPartitionBase::createLocalDofOrderings: " << dofNosLocal_;

  // create IS (indexSet) dofNosLocalIS_;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(PETSC_COMM_SELF, nDofsLocalWithGhosts, dofNosLocal_.data(), PETSC_COPY_VALUES, &dofNosLocalIS_); CHKERRV(ierr);

  // create IS dofNosLocalNonGhostIS_;
  ierr = ISCreateGeneral(PETSC_COMM_SELF, nDofsLocalWithoutGhosts, dofNosLocal_.data(), PETSC_COPY_VALUES, &dofNosLocalNonGhostIS_); CHKERRV(ierr);
}

const std::vector<PetscInt> &MeshPartitionBase::dofNosLocal() const
{
  return dofNosLocal_;
}

const std::vector<dof_no_t> &MeshPartitionBase::dofNosLocalNaturalOrdering() const
{
  return dofNosLocal_;
}

//! get a PETSc IS (index set) with the same information as dofNosLocal_
const IS &MeshPartitionBase::dofNosLocalIS() const
{
  return dofNosLocalIS_;
}

const IS &MeshPartitionBase::dofNosLocalNonGhostIS() const
{
  return dofNosLocalNonGhostIS_;
}

}  // namespace
