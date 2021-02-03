#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

#include "utility/mpi_utility.h"

// ---- incompressible case ----
//! constructor
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
PartitionedPetscVecForHyperelasticity(
  std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>> meshPartitionDisplacements,
  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure,
  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpaceType,nComponents>> dirichletBoundaryConditions, std::string name) :
  PartitionedPetscVecWithDirichletBc<DisplacementsFunctionSpaceType,nComponents+1,nComponents>(meshPartitionDisplacements, dirichletBoundaryConditions, name, true),
  meshPartitionPressure_(meshPartitionPressure)
{
  componentNoPressure_ = nComponents;   // the pressure component is the last one, equal to nComponents (=3 or 6)

  // initialize variables for 3 or 6 displacement/velocity components and 1 pressure component
  LOG(DEBUG) << "\"" << this->name_ << "\" initialize PartitionedPetscVecForHyperelasticity";
  this->initialize(meshPartitionPressure_->nDofsLocalWithoutGhosts());

  // initialize last component for pressure
  initializeForPressure();

  // create the Petsc Vec
  this->createVector();
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
void PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
initializeForPressure()
{
  VLOG(1) << "initializeForPressure, nEntriesGlobal_: " << this->nEntriesGlobal_ << ", nEntriesLocal_: " << this->nEntriesLocal_
   << ", nNonBcDofsGhosts_: " << this->nNonBcDofsGhosts_;

  const int nDofsLocalWithGhosts = meshPartitionPressure_->nDofsLocalWithGhosts();
  const int nDofsLocalWithoutGhosts = meshPartitionPressure_->nDofsLocalWithoutGhosts();

  VLOG(1) << "initializeForPressure nDofsLocalWithGhosts: " << nDofsLocalWithGhosts << ", nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts;

  // the local without ghosts number of entries in the vector, without the Dirichlet BC dofs
  this->nNonBcDofsWithoutGhosts_[componentNoPressure_] = nDofsLocalWithoutGhosts;

  // number of ghost values
  int nGhostsPressure = nDofsLocalWithGhosts - nDofsLocalWithoutGhosts;

  this->nonBcGhostDofNosGlobal_.reserve(this->nonBcGhostDofNosGlobal_.size() + nGhostsPressure);

  // set mapping from component no and local dof no to the global numbering, for non-ghost dofs
  this->dofNoLocalToDofNoNonBcGlobal_[componentNoPressure_].resize(nDofsLocalWithoutGhosts);
  std::iota(this->dofNoLocalToDofNoNonBcGlobal_[componentNoPressure_].begin(), this->dofNoLocalToDofNoNonBcGlobal_[componentNoPressure_].end(),
            this->nonBcDofNoGlobalBegin_+this->nDofsLocal_);
  this->dofNoLocalToDofNoNonBcGlobal_[componentNoPressure_].reserve(nDofsLocalWithGhosts);

  // add ghost numbers for pressure
  const std::vector<PetscInt> &ghostDofNosGlobalPetsc = meshPartitionPressure_->ghostDofNosGlobalPetsc();

  int nRanks = this->meshPartitionPressure_->nRanks();
  std::vector<int> nDisplacementDofs(nRanks);

  // distribute the sizes on the processes
  int nPreviousDisplacementDofs = 0;
  MPIUtility::handleReturnValue(MPI_Scan(&this->nDofsLocal_, &nPreviousDisplacementDofs, 1, MPI_INT, MPI_SUM, this->meshPartitionPressure_->mpiCommunicator()), "MPI_Scan");
  MPIUtility::handleReturnValue(MPI_Allgather(&nPreviousDisplacementDofs, 1, MPI_INT, nDisplacementDofs.data(), 1, MPI_INT, this->meshPartitionPressure_->mpiCommunicator()), "MPI_Allgather");

  VLOG(1) << "nPreviousDisplacementDofs: " << nPreviousDisplacementDofs << ", nDisplacementDofs: " << nDisplacementDofs;

  int i = 0;
  for (int dofNoLocal = nDofsLocalWithoutGhosts; dofNoLocal < nDofsLocalWithGhosts; dofNoLocal++, i++)
  {
    node_no_t nodeNoLocal = dofNoLocal / PressureFunctionSpaceType::nDofsPerNode();

    int neighbourRankNo;
    meshPartitionPressure_->isNonGhost(nodeNoLocal, neighbourRankNo);

    global_no_t nonBcGhostDofNoGlobal = ghostDofNosGlobalPetsc[i] + nDisplacementDofs[neighbourRankNo];
    this->nonBcGhostDofNosGlobal_.push_back(nonBcGhostDofNoGlobal);

    this->dofNoLocalToDofNoNonBcGlobal_[componentNoPressure_].push_back(nonBcGhostDofNoGlobal);

    VLOG(1) << "dofNoLocal " << dofNoLocal << " is on rank " << neighbourRankNo << ", this rank has previous nDisplacementDofs: "
      << nDisplacementDofs[neighbourRankNo]
      << ", ghost dof globalPetsc: " << ghostDofNosGlobalPetsc[i] << ", nonBcGhostDofNoGlobal: " << nonBcGhostDofNoGlobal;
  }

  // set mapping from component no and local dof no to the local number of the non-bc dof numbering
  this->dofNoLocalToDofNoNonBcLocal_[componentNoPressure_].resize(nDofsLocalWithGhosts);

  // non-ghost dofs
  std::iota(this->dofNoLocalToDofNoNonBcLocal_[componentNoPressure_].begin(), this->dofNoLocalToDofNoNonBcLocal_[componentNoPressure_].begin()+nDofsLocalWithoutGhosts,
            this->nEntriesLocal_-nDofsLocalWithoutGhosts);

  // ghost dofs
  std::iota(this->dofNoLocalToDofNoNonBcLocal_[componentNoPressure_].begin()+nDofsLocalWithoutGhosts, this->dofNoLocalToDofNoNonBcLocal_[componentNoPressure_].end(),
            this->nEntriesLocal_+this->nNonBcDofsGhosts_);

  this->nNonBcDofsGhosts_ += nGhostsPressure;

  // there are no boundary conditions for pressure component
  this->boundaryConditionValues_[componentNoPressure_].resize(nDofsLocalWithGhosts, -1.0);
  this->isPrescribed_[componentNoPressure_].resize(nDofsLocalWithGhosts, false);

  VLOG(1) << "after initializeForPressure:";
  VLOG(1) << "dofNoLocalToDofNoNonBcLocal_: " << this->dofNoLocalToDofNoNonBcLocal_;
  VLOG(1) << "dofNoLocalToDofNoNonBcGlobal_: " << this->dofNoLocalToDofNoNonBcGlobal_;
  VLOG(1) << "nonBcGhostDofNosGlobal_: " << this->nonBcGhostDofNosGlobal_;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
meshPartitionPressure()
{
  return meshPartitionPressure_;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
void PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
dumpGlobalNatural(std::string filename)
{
  // write file
  std::ofstream file;
  std::string vectorName = filename;
  if (vectorName.find("/") != std::string::npos)
    vectorName = vectorName.substr(vectorName.rfind("/")+1);
  filename += std::string(".m");
  OutputWriter::Generic::openFile(file, filename);

  // write header
  file << "% " << this->name_ << ", " << this->nEntriesGlobal_ << " entries, " << this->meshPartition_->nRanks() << " MPI ranks" << std::endl;
  file << getString(false, vectorName);
  LOG(INFO) << "Vector \"" << filename << "\" written.";
  file.close();
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
displacementDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of displacement dofs
  int nDisplacementDofsLocal = nDisplacementDofsWithoutBcLocal();
  std::vector<dof_no_t> indices(nDisplacementDofsLocal);
  std::iota(indices.begin(), indices.end(), this->nonBcDofNoGlobalBegin_);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nDisplacementDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
velocityDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of displacement dofs
  int nDisplacementDofsLocal = nDisplacementDofsWithoutBcLocal();
  int nVelocityDofsLocal = nVelocityDofsWithoutBcLocal();
  std::vector<dof_no_t> indices(nVelocityDofsLocal);
  std::iota(indices.begin(), indices.end(), this->nonBcDofNoGlobalBegin_ + nDisplacementDofsLocal);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nVelocityDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
pressureDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of pressure dofs
  int nPressureDofsLocal = this->nNonBcDofsWithoutGhosts_[componentNoPressure_];
  std::vector<dof_no_t> indices(nPressureDofsLocal);
  int startIndex = this->nonBcDofNoGlobalBegin_ + nDisplacementDofsWithoutBcLocal() + nVelocityDofsWithoutBcLocal();
  std::iota(indices.begin(), indices.end(), startIndex);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nPressureDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
dof_no_t PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
nDisplacementDofsWithoutBcLocal()
{
  dof_no_t result = 0;
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    result += this->nNonBcDofsWithoutGhosts_[componentNo];
  }
  return result;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
dof_no_t PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
nVelocityDofsWithoutBcLocal()
{
  dof_no_t result = 0;
  for (int componentNo = 3; componentNo < nComponents; componentNo++)
  {
    result += this->nNonBcDofsWithoutGhosts_[componentNo];
  }
  return result;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
bool PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>::
containsNanOrInf()
{
  // get all local values
  std::vector<double> values;

  for (int componentNo = 0; componentNo < nComponents-1; componentNo++)
  {
    values.resize(this->meshPartition_->nDofsLocalWithoutGhosts());

    // void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const;
    this->getValues(componentNo, this->meshPartition_->nDofsLocalWithoutGhosts(), this->meshPartition_->dofNosLocal().data(), values.data());

    // loop over values and check if they are neither nan nor inf
    for (int i = 0; i < values.size(); i++)
    {
      if (!std::isfinite(values[i]) || fabs(values[i]) > 1e+75)
      {
        LOG(ERROR) << "containsNanOrInf(): value " << i << "/" << values.size() << ", component " << componentNo << "/" << nComponents << ": " << values[i];
        return true;
      }
    }
  }
  return false;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term> &vector)
{
  vector.output(stream);
  stream << vector.getString();
  return stream;
}
