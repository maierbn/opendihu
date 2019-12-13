#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

#include "utility/mpi_utility.h"

//! constructor
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
PartitionedPetscVecForHyperelasticity(
  std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>> meshPartitionDisplacements,
  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure,
  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpaceType,nComponents>> dirichletBoundaryConditions, std::string name) :
  PartitionedPetscVecWithDirichletBc<DisplacementsFunctionSpaceType,nComponents+1,nComponents>(meshPartitionDisplacements, dirichletBoundaryConditions, name), meshPartitionPressure_(meshPartitionPressure)
{
  componentNoPressure_ = nComponents;   // the pressure component is the last one, equal to nComponents (=3 or 6)

  // initialize variables for 3 or 6 displacement/velocity components
  this->initialize(nComponents, meshPartitionPressure_->nDofsLocalWithoutGhosts());

  // initialize last component for pressure
  initializeForPressure();

  // create the Petsc Vec
  this->createVector();

  // debugging output
  /*std::stringstream stream;
  this->output(stream);
  LOG(DEBUG) << stream.str();*/
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
void PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
initializeForPressure()
{
  VLOG(1) << "initializeForPressure, nEntriesGlobal_: " << this->nEntriesGlobal_ << ", nEntriesLocal_: " << this->nEntriesLocal_;

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
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
meshPartitionPressure()
{
  return meshPartitionPressure_;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
void PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
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

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
std::string PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
getString(bool horizontal, std::string vectorName) const
{
  // do not assemble a horizontal string for console in release mode, because this is only needed for debugging output
#ifdef NDEBUG
  if (horizontal)
    return std::string("");
#endif

  std::stringstream result;
  int ownRankNo = this->meshPartition_->ownRankNo();

  // the master rank collects all data and writes the file
  if (ownRankNo == 0)
  {
    // handle displacement values
    std::vector<std::vector<global_no_t>> dofNosGlobalNatural(this->meshPartition_->nRanks());
    std::vector<std::vector<double>> values(this->meshPartition_->nRanks());
    std::vector<int> nDofsOnRank(this->meshPartition_->nRanks());
    std::vector<int> nPressureDofsOnRank(this->meshPartition_->nRanks());

    VLOG(1) << "values: " << values;

    // handle own rank
    // get dof nos in global natural ordering
    int nDofsLocalWithoutGhosts = this->meshPartition_->nDofsLocalWithoutGhosts();
    nDofsOnRank[0] = nDofsLocalWithoutGhosts;

    std::vector<global_no_t> dofNosGlobalNaturalOwn;
    this->meshPartition_->getDofNosGlobalNatural(dofNosGlobalNaturalOwn);

    dofNosGlobalNatural[0].resize(nDofsLocalWithoutGhosts);
    std::copy(dofNosGlobalNaturalOwn.begin(), dofNosGlobalNaturalOwn.end(), dofNosGlobalNatural[0].begin());

    // get displacement values
    values[0].resize(nComponents*nDofsLocalWithoutGhosts);
    for (int i = 0; i < nComponents; i++)
    {
      this->getValues(i, nDofsLocalWithoutGhosts, this->meshPartition_->dofNosLocal().data(), values[0].data() + i*nDofsLocalWithoutGhosts);
    }

    VLOG(1) << "get own displacement values: " << values;

    for (int rankNoK = 0; rankNoK < this->meshPartition_->nRanks(2); rankNoK++)
    {
      for (int rankNoJ = 0; rankNoJ < this->meshPartition_->nRanks(1); rankNoJ++)
      {
        for (int rankNoI = 0; rankNoI < this->meshPartition_->nRanks(0); rankNoI++)
        {
          int rankNo = rankNoK*this->meshPartition_->nRanks(1)*this->meshPartition_->nRanks(0)
            + rankNoJ*this->meshPartition_->nRanks(0) + rankNoI;

          VLOG(1) << "rank " << rankNo;

          // determine number of dofs on rank
          nDofsOnRank[rankNo] = this->meshPartition_->nNodesLocalWithoutGhosts(2, rankNoK)
          * this->meshPartition_->nNodesLocalWithoutGhosts(1, rankNoJ)
          * this->meshPartition_->nNodesLocalWithoutGhosts(0, rankNoI)
          * DisplacementsFunctionSpaceType::nDofsPerNode();

          if (rankNo == 0)
            continue;

          VLOG(1) << ", nDofsOnRank: " << nDofsOnRank[rankNo];

          VLOG(1) << "recv from " << rankNo << " " << nDofsOnRank[rankNo] << " dofs and displacements values";

          // receive dof nos
          dofNosGlobalNatural[rankNo].resize(nDofsOnRank[rankNo]);
          MPIUtility::handleReturnValue(MPI_Recv(dofNosGlobalNatural[rankNo].data(), nDofsOnRank[rankNo], MPI_UNSIGNED_LONG_LONG,
                                                rankNo, 0, this->meshPartition_->rankSubset()->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          VLOG(1) << "received displacements dofs: " << dofNosGlobalNatural[rankNo];

          // receive values
          values[rankNo].resize(nComponents*nDofsOnRank[rankNo]);
          MPIUtility::handleReturnValue(MPI_Recv(values[rankNo].data(), nComponents*nDofsOnRank[rankNo], MPI_DOUBLE,
                                                rankNo, 0, this->meshPartition_->rankSubset()->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          VLOG(1) << "received displacements values: " << values[rankNo];
        }
      }
    }

    VLOG(1) << "values: " << values;

    // handle pressure values
    std::vector<std::vector<global_no_t>> dofNosGlobalNaturalPressure(this->meshPartition_->nRanks());
    std::vector<std::vector<double>> valuesPressure(this->meshPartition_->nRanks());

    // handle own rank
    // get global natural dof nos for pressure
    nDofsLocalWithoutGhosts = this->meshPartitionPressure_->nDofsLocalWithoutGhosts();
    nPressureDofsOnRank[0] = nDofsLocalWithoutGhosts;

    VLOG(1) << "nRanks: " << this->meshPartition_->nRanks() << ", nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts;

    std::vector<global_no_t> dofNosGlobalNaturalPressureOwn;
    this->meshPartitionPressure_->getDofNosGlobalNatural(dofNosGlobalNaturalPressureOwn);

    dofNosGlobalNaturalPressure[0].resize(nDofsLocalWithoutGhosts);
    std::copy(dofNosGlobalNaturalPressureOwn.begin(), dofNosGlobalNaturalPressureOwn.end(), dofNosGlobalNaturalPressure[0].begin());

    // get pressure values
    valuesPressure[0].resize(nDofsLocalWithoutGhosts);
    this->getValues(componentNoPressure_, nDofsLocalWithoutGhosts, this->meshPartitionPressure_->dofNosLocal().data(), valuesPressure[0].data());

    // loop over other ranks
    for (int rankNoK = 0; rankNoK < this->meshPartitionPressure_->nRanks(2); rankNoK++)
    {
      for (int rankNoJ = 0; rankNoJ < this->meshPartitionPressure_->nRanks(1); rankNoJ++)
      {
        for (int rankNoI = 0; rankNoI < this->meshPartitionPressure_->nRanks(0); rankNoI++)
        {
          int rankNo = rankNoK*this->meshPartitionPressure_->nRanks(1)*this->meshPartitionPressure_->nRanks(0)
            + rankNoJ*this->meshPartitionPressure_->nRanks(0) + rankNoI;

          // determine number of dofs on rank
          nPressureDofsOnRank[rankNo] = this->meshPartitionPressure_->nNodesLocalWithoutGhosts(2, rankNoK)
          * this->meshPartitionPressure_->nNodesLocalWithoutGhosts(1, rankNoJ)
          * this->meshPartitionPressure_->nNodesLocalWithoutGhosts(0, rankNoI)
          * DisplacementsFunctionSpaceType::nDofsPerNode();

          VLOG(1) << "pressure rank " << rankNo << " n dofs: " << nPressureDofsOnRank[rankNo];

          if (rankNo == 0)
            continue;

          VLOG(1) << "recv from " << rankNo << " " << nPressureDofsOnRank[rankNo] << " dofs and pressure values";

          // receive dof nos
          dofNosGlobalNaturalPressure[rankNo].resize(nPressureDofsOnRank[rankNo]);
          MPIUtility::handleReturnValue(MPI_Recv(dofNosGlobalNaturalPressure[rankNo].data(), nPressureDofsOnRank[rankNo], MPI_UNSIGNED_LONG_LONG,
                                                rankNo, 0, this->meshPartitionPressure_->rankSubset()->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          VLOG(1) << "received pressureDofs: " << dofNosGlobalNaturalPressure[rankNo];

          // receive values
          valuesPressure[rankNo].resize(nPressureDofsOnRank[rankNo]);
          MPIUtility::handleReturnValue(MPI_Recv(valuesPressure[rankNo].data(), nPressureDofsOnRank[rankNo], MPI_DOUBLE,
                                                rankNo, 0, this->meshPartitionPressure_->rankSubset()->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          VLOG(1) << "received pressure values: " << valuesPressure[rankNo];
        }
      }
    }

    // sort displacement values according to global natural dof no
    std::vector<std::pair<global_no_t,std::array<double,nComponents>>> displacementEntries;
    displacementEntries.reserve(this->nEntriesGlobal_);

    for (int rankNo = 0; rankNo < this->meshPartition_->nRanks(); rankNo++)
    {
      assert(nDofsOnRank[rankNo] == dofNosGlobalNatural[rankNo].size());
      for (int i = 0; i < nDofsOnRank[rankNo]; i++)
      {
        std::pair<global_no_t,std::array<double,nComponents>> displacementEntry;
        displacementEntry.first = dofNosGlobalNatural[rankNo][i];

        std::array<double,nComponents> valuesDof;
        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          valuesDof[componentNo] = values[rankNo][componentNo*nDofsOnRank[rankNo] + i];
        }
        displacementEntry.second = valuesDof;

        displacementEntries.push_back(displacementEntry);
      }
    }

    // sort list according to dof no.s
    std::sort(displacementEntries.begin(), displacementEntries.end(),
              [&](std::pair<global_no_t,std::array<double,nComponents>> a, std::pair<global_no_t,std::array<double,nComponents>> b)
    {
      return a.first < b.first;
    });

    // sort pressure values according to global natural dof no
    std::vector<std::pair<global_no_t,double>> pressureEntries;
    pressureEntries.reserve(this->nEntriesGlobal_);

    for (int rankNo = 0; rankNo < this->meshPartitionPressure_->nRanks(); rankNo++)
    {
      int nDofsOnRank = dofNosGlobalNaturalPressure[rankNo].size();
      for (int i = 0; i < nDofsOnRank; i++)
      {
        pressureEntries.push_back(std::pair<global_no_t,double>(
          dofNosGlobalNaturalPressure[rankNo][i],
          valuesPressure[rankNo][i]
        ));
      }
    }

    // sort list according to dof no.s
    std::sort(pressureEntries.begin(), pressureEntries.end(),
              [&](std::pair<global_no_t,double> a, std::pair<global_no_t,double> b)
    {
      return a.first < b.first;
    });


    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "dofNosGlobalNatural: " << dofNosGlobalNatural;
      VLOG(1) << "values: " << values;
      VLOG(1) << "nDofsOnRank: " << nDofsOnRank;
      VLOG(1) << "nPressureDofsOnRank: " << nPressureDofsOnRank;
      VLOG(1) << "valuesPressure: " << valuesPressure;
      VLOG(1) << "displacementEntries: " << displacementEntries;
      VLOG(1) << "pressureEntries: " << pressureEntries;
    }


    // write file
    std::string newline = "\n";
    std::string separator = ", ";
    std::array<std::string,nComponents> componentNames;

    componentNames[0] = "ux";
    componentNames[1] = "uy";
    componentNames[2] = "uz";
    if (nComponents == 6)
    {
      componentNames[3] = "vx";
      componentNames[4] = "vy";
      componentNames[5] = "vz";
    }

    if (horizontal)
    {
      newline = "";
      separator = ", ";
      result << std::endl;
    }
    else
    {
      result << vectorName << "r" << this->meshPartitionPressure_->nRanks() << " = [";
    }

    // loop over not-pressure components (u and possibly v)
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      // start of component
      if (horizontal)
      {
        result << componentNames[componentNo] << " = [" << newline;  // print e.g. "ux = ["
      }

      // write displacement values
      for (int i = 0; i < displacementEntries.size(); i++)
      {
        if (i != 0)
          result << separator;
        result << displacementEntries[i].second[0];
      }

      // end of component
      if (horizontal)
      {
        result << newline << "]; " << std::endl;
      }
      else
      {
        result << ", ...\n ";
      }
    }

    if (horizontal)
    {
      result << " p = [" << newline;
    }
    else
    {
      result << ", ...\n ";
    }

    for (int i = 0; i < pressureEntries.size(); i++)
    {
      if (i != 0)
        result << separator;
      result << pressureEntries[i].second;
    }

    if (horizontal)
    {
      result << newline << "];" << std::endl;
    }
    else
    {
      result << "]; " << std::endl;
    }
  }
  else
  {
    // send global natural dof nos for displacements
    int nDofsLocalWithoutGhosts = this->meshPartition_->nDofsLocalWithoutGhosts();

    std::vector<global_no_t> dofNosGlobalNatural;
    this->meshPartition_->getDofNosGlobalNatural(dofNosGlobalNatural);

    assert(dofNosGlobalNatural.size() == nDofsLocalWithoutGhosts);

    VLOG(1) << "send to 0 " << nDofsLocalWithoutGhosts << " dofs and displacements values";

    MPIUtility::handleReturnValue(MPI_Send(dofNosGlobalNatural.data(), nDofsLocalWithoutGhosts, MPI_UNSIGNED_LONG_LONG,
                                           0, 0, this->meshPartition_->rankSubset()->mpiCommunicator()), "MPI_Send");

    VLOG(1) << "sent displacements dofs: " << dofNosGlobalNatural;

    // send displacement values
    std::vector<double> values(nComponents*nDofsLocalWithoutGhosts);

    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      this->getValues(componentNo, nDofsLocalWithoutGhosts, this->meshPartition_->dofNosLocal().data(), values.data() + componentNo*nDofsLocalWithoutGhosts);
    }

    MPIUtility::handleReturnValue(MPI_Send(values.data(), nComponents*nDofsLocalWithoutGhosts, MPI_DOUBLE,
                                           0, 0, this->meshPartition_->rankSubset()->mpiCommunicator()), "MPI_Send");

    VLOG(1) << "sent displacements values: " << values;

    // send global natural dof nos for pressure
    nDofsLocalWithoutGhosts = this->meshPartitionPressure_->nDofsLocalWithoutGhosts();

    VLOG(1) << "send to 0 " << nDofsLocalWithoutGhosts << " pressure dofs and values";

    std::vector<global_no_t> dofNosGlobalNaturalPressure;
    this->meshPartitionPressure_->getDofNosGlobalNatural(dofNosGlobalNaturalPressure);

    assert(dofNosGlobalNaturalPressure.size() == nDofsLocalWithoutGhosts);

    MPIUtility::handleReturnValue(MPI_Send(dofNosGlobalNaturalPressure.data(), nDofsLocalWithoutGhosts, MPI_UNSIGNED_LONG_LONG,
                                           0, 0, this->meshPartitionPressure_->rankSubset()->mpiCommunicator()), "MPI_Send");

    VLOG(1) << "sent pressure dofs: " << dofNosGlobalNaturalPressure;

    // send pressure values
    std::vector<double> valuesPressure(nDofsLocalWithoutGhosts);

    this->getValues(componentNoPressure_, nDofsLocalWithoutGhosts, this->meshPartitionPressure_->dofNosLocal().data(), valuesPressure.data());

    MPIUtility::handleReturnValue(MPI_Send(valuesPressure.data(), nDofsLocalWithoutGhosts, MPI_DOUBLE,
                                           0, 0, this->meshPartitionPressure_->rankSubset()->mpiCommunicator()), "MPI_Send");

    VLOG(1) << "sent pressure values: " << valuesPressure;
  }
  return result.str();
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
displacementDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of displacement dofs
  int nDisplacementDofsLocal = nDisplacementDofsWithoutBcLocal();
  std::vector<int> indices(nDisplacementDofsLocal);
  std::iota(indices.begin(), indices.end(), this->nonBcDofNoGlobalBegin_);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nDisplacementDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
velocityDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of displacement dofs
  int nDisplacementDofsLocal = nDisplacementDofsWithoutBcLocal();
  int nVelocityDofsLocal = nVelocityDofsWithoutBcLocal();
  std::vector<int> indices(nVelocityDofsLocal);
  std::iota(indices.begin(), indices.end(), this->nonBcDofNoGlobalBegin_ + nDisplacementDofsLocal);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nVelocityDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
pressureDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of pressure dofs
  int nPressureDofsLocal = this->nNonBcDofsWithoutGhosts_[componentNoPressure_];
  std::vector<int> indices(nPressureDofsLocal);
  int startIndex = this->nonBcDofNoGlobalBegin_ + nDisplacementDofsWithoutBcLocal() + nVelocityDofsWithoutBcLocal();
  std::iota(indices.begin(), indices.end(), startIndex);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nPressureDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
dof_no_t PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
nDisplacementDofsWithoutBcLocal()
{
  dof_no_t result = 0;
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    result += this->nNonBcDofsWithoutGhosts_[componentNo];
  }
  return result;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
dof_no_t PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nComponents>::
nVelocityDofsWithoutBcLocal()
{
  dof_no_t result = 0;
  for (int componentNo = 3; componentNo < nComponents; componentNo++)
  {
    result += this->nNonBcDofsWithoutGhosts_[componentNo];
  }
  return result;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType> &vector)
{
  vector.output(stream);
  stream << vector.getString();
  return stream;
}
