#include "partition/00_mesh_partition_base.h"

#include <algorithm>

#include "utility/string_utility.h"
#include "semt/Semt.h"
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

MPI_Comm MeshPartitionBase::mpiCommunicator() const
{
  return rankSubset_->mpiCommunicator();
}

  //! fill the dofLocalNo vectors
void MeshPartitionBase::createLocalDofOrderings(dof_no_t nDofsLocal)
{
  // fill dofNosLocal_ 
  dofNosLocal_.resize(nDofsLocal);
  std::iota(dofNosLocal_.begin(), dofNosLocal_.end(), 0);
  
  
  LOG(DEBUG) << " MeshPartitionBase::createLocalDofOrderings: " << dofNosLocal_;
  
  // create IS (indexSet) dofNosLocalIS_;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(PETSC_COMM_SELF, dofNosLocal_.size(), dofNosLocal_.data(), PETSC_COPY_VALUES, &dofNosLocalIS_); CHKERRV(ierr);
}

const std::vector<PetscInt> &MeshPartitionBase::dofNosLocal() const
{
  return dofNosLocal_;
}

//! get a PETSc IS (index set) with the same information as dofNosLocal_
const IS &MeshPartitionBase::dofNosLocalIS() const
{
  return dofNosLocalIS_;
}

}  // namespace
