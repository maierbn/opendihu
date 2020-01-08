#include "partition/partitioned_petsc_vec/01_partitioned_petsc_vec_with_dirichlet_bc.h"

#include "utility/mpi_utility.h"

//! constructor
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
PartitionedPetscVecWithDirichletBc(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition,
                                   std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,nComponentsDirichletBc>> dirichletBoundaryConditions, std::string name) :
  PartitionedPetscVec<FunctionSpaceType,nComponents>(meshPartition, name), dirichletBoundaryConditions_(dirichletBoundaryConditions)
{
  initialize();
  createVector();
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
global_no_t PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nonBCDofNoGlobal(int componentNo, dof_no_t localDofNo) const
{
  return dofNoLocalToDofNoNonBcGlobal_[componentNo][localDofNo];
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
dof_no_t PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nonBCDofNoLocal(int componentNo, dof_no_t localDofNo) const
{
  return dofNoLocalToDofNoNonBcLocal_[componentNo][localDofNo];
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
bool PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
isPrescribed(int componentNo, dof_no_t localDofNo) const
{
  return isPrescribed_[componentNo][localDofNo];
}

//! get a reference to the internal dofNoLocalToDofNoNonBcGlobal_ data structure
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
const std::array<std::vector<dof_no_t>,nComponents> &PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
dofNoLocalToDofNoNonBcGlobal() const
{
  return dofNoLocalToDofNoNonBcGlobal_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
initialize(int offsetInGlobalNumberingPerRank)
{
  LOG(DEBUG) << "\"" << this->name_ << "\" PartitionedPetscVecWithDirichletBc createVector with " << nComponentsDirichletBc << " components, size local: " << this->meshPartition_->nNodesLocalWithoutGhosts()
    << ", global: " << this->meshPartition_->nNodesGlobal() << ", offsetInGlobalNumberingPerRank: " << offsetInGlobalNumberingPerRank
    << ", ghost dof nos global/petsc: " << this->meshPartition_->ghostDofNosGlobalPetsc();

  // normally, nComponentsDirichletBc == nComponents
  // only the deriving class can set nComponents > nComponentsDirichletBc, then it has to also properly initialize everything above nComponentsDirichletBc

  // define abbreviation variables
  const int nDofsLocalWithoutGhosts = this->meshPartition_->nDofsLocalWithoutGhosts();
  const int nDofsLocalWithGhosts = this->meshPartition_->nDofsLocalWithGhosts();

  // determine the number of local entries of the vector
  nEntriesLocal_ = 0;

  // loop over components
  for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
  {
    int nBoundaryConditionDofs = 0;
    nBoundaryConditionDofs = dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].dofNosLocal.size();
    nNonBcDofsWithoutGhosts_[componentNo] = nDofsLocalWithoutGhosts - nBoundaryConditionDofs;
    nEntriesLocal_ += nNonBcDofsWithoutGhosts_[componentNo];
  }

  // add space for pressure components, this is only !=0 for PartitionedPetscVecForHyperelasticity
  nDofsLocal_ = nEntriesLocal_;
  nEntriesLocal_ += offsetInGlobalNumberingPerRank;

  // determine global number of non-BC dofs over all components
  nEntriesGlobal_ = 0;
  MPIUtility::handleReturnValue(MPI_Allreduce(&nEntriesLocal_, &nEntriesGlobal_, 1, MPI_INT, MPI_SUM, this->meshPartition_->mpiCommunicator()), "MPI_Allreduce");

  // get number of non-bc dofs on previous ranks
  nonBcDofNoGlobalBegin_ = 0;
  MPIUtility::handleReturnValue(MPI_Exscan(&nEntriesLocal_, &nonBcDofNoGlobalBegin_, 1, MPI_INT, MPI_SUM, this->meshPartition_->mpiCommunicator()), "MPI_Exscan");

  // setup non-bc numberings (local and global) for non-ghost dofs, also store boundaryConditionValues_ and isPrescribed_ for non-ghost dofs
  global_no_t dofNoNonBcGlobal = nonBcDofNoGlobalBegin_;

  // loop over components
  for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
  {
    std::vector<dof_no_t>::const_iterator boundaryConditionDofLocalNosIter = dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].dofNosLocal.begin();
    std::vector<double>::const_iterator boundaryConditionValuesIter = dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].values.begin();

    dofNoLocalToDofNoNonBcGlobal_[componentNo].resize(nDofsLocalWithGhosts);
    dofNoLocalToDofNoNonBcLocal_[componentNo].resize(nDofsLocalWithGhosts);
    boundaryConditionValues_[componentNo].resize(nDofsLocalWithGhosts);
    isPrescribed_[componentNo].resize(nDofsLocalWithGhosts, false);

    // loop over local dofs
    for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
    {
      VLOG(1) << "componentNo " << componentNo << " dofǸoLocal " << dofNoLocal;

      // if dof has bc
      if (boundaryConditionDofLocalNosIter != dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].dofNosLocal.end())
      {
        VLOG(1) << " next bc: " << *boundaryConditionDofLocalNosIter;

        if (*boundaryConditionDofLocalNosIter == dofNoLocal)
        {
          VLOG(1) << " -> is bc";

          dofNoLocalToDofNoNonBcGlobal_[componentNo][dofNoLocal] = -1;
          dofNoLocalToDofNoNonBcLocal_[componentNo][dofNoLocal] = -1;
          boundaryConditionValues_[componentNo][dofNoLocal] = *boundaryConditionValuesIter;
          isPrescribed_[componentNo][dofNoLocal] = true;

          boundaryConditionDofLocalNosIter++;
          boundaryConditionValuesIter++;
          continue;
        }
      }
      VLOG(1) << " is no bc, assign new no. " << dofNoNonBcGlobal;

      // here, dofNoLocal is a non-BC dof
      dofNoLocalToDofNoNonBcGlobal_[componentNo][dofNoLocal] = dofNoNonBcGlobal;
      dofNoLocalToDofNoNonBcLocal_[componentNo][dofNoLocal] = dofNoNonBcGlobal - nonBcDofNoGlobalBegin_;
      boundaryConditionValues_[componentNo][dofNoLocal] = -1.0;

      dofNoNonBcGlobal++;
    }
  }

  VLOG(1) << "dofNoLocalToDofNoNonBcGlobal_ (only non-ghost dofs are set so far): " << dofNoLocalToDofNoNonBcGlobal_;
  VLOG(1) << "dofNoLocalToDofNoNonBcLocal_ (only non-ghost dofs are set so far): " << dofNoLocalToDofNoNonBcLocal_;

  // next, we communicate with neighbors to obtain non-bc numbering for ghost dofs

  // get vector of ghost dofs that are needed from neighbouring ranks
  struct DofsRequest
  {
    std::vector<global_no_t> dofNosGlobalPetsc;
    std::vector<dof_no_t> dofNosLocal;
  };

  std::map<int, DofsRequest> requestDofsFromRanks;

  for (dof_no_t dofNoLocal = nDofsLocalWithoutGhosts; dofNoLocal < this->meshPartition_->nDofsLocalWithGhosts(); dofNoLocal++)
  {
    global_no_t dofNoGlobalPetsc = this->meshPartition_->getDofNoGlobalPetsc(dofNoLocal);

    // if the current dof is a ghost dof on the own rank (should be)
    node_no_t nodeNoLocal = dofNoLocal / FunctionSpaceType::nDofsPerNode();
    int neighbourRankNo = 0;
    if (!this->meshPartition_->isNonGhost(nodeNoLocal, neighbourRankNo))
    {
      requestDofsFromRanks[neighbourRankNo].dofNosLocal.push_back(dofNoLocal);
      requestDofsFromRanks[neighbourRankNo].dofNosGlobalPetsc.push_back(dofNoGlobalPetsc);
    }
    else
    {
      LOG(FATAL) << "ghost dof not recognized as ghost dof (isNonGhost is errorneous)";
    }
  }

  // exchange, how many values should be sent to which rank
  VLOG(1) << "rankSubset " << *this->meshPartition_->rankSubset() << ", create new window";

  int nRanks = this->meshPartition_->nRanks();
  int ownRankNo = this->meshPartition_->ownRankNo();

  // create remote accessible memory
  std::vector<int> remoteAccessibleMemory(nRanks, 0);
  int nBytes = nRanks * sizeof(int);
  int displacementUnit = sizeof(int);
  MPI_Win mpiMemoryWindow;
  MPIUtility::handleReturnValue(MPI_Win_create((void *)remoteAccessibleMemory.data(), nBytes, displacementUnit, MPI_INFO_NULL, this->meshPartition_->mpiCommunicator(), &mpiMemoryWindow), "MPI_Win_create");

  std::vector<int> localMemory(nRanks);

  // put number of requested ghost dofs to the corresponding processes
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks)
  {
    int foreignRankNo = requestDofsFromRank.first;
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();
    localMemory[foreignRankNo] = nRequestedDofs;

    VLOG(1) << "put value " << localMemory[foreignRankNo] << " to rank " << foreignRankNo << ", offset " << ownRankNo;

    // start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
    MPIUtility::handleReturnValue(MPI_Win_lock(MPI_LOCK_SHARED, foreignRankNo, 0, mpiMemoryWindow), "MPI_Win_lock");

    MPIUtility::handleReturnValue(MPI_Put(&localMemory[foreignRankNo], 1, MPI_INT, foreignRankNo, ownRankNo, 1, MPI_INT, mpiMemoryWindow), "MPI_Put");

    MPIUtility::handleReturnValue(MPI_Win_unlock(foreignRankNo, mpiMemoryWindow), "MPI_Win_unlock");

  }

  MPIUtility::handleReturnValue(MPI_Win_fence(MPI_MODE_NOSUCCEED, mpiMemoryWindow), "MPI_Win_fence");

  std::vector<std::pair<int,int>> nDofRequestedFromRanks;   /// (foreignRank,nDofs), number of dofs requested by and to be send to foreignRank
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    VLOG(1) << " rank " << rankNo << " nRequestedDofs: " << remoteAccessibleMemory[rankNo];
    if (remoteAccessibleMemory[rankNo] > 0)
    {
      nDofRequestedFromRanks.push_back(std::pair<int,int>(rankNo,remoteAccessibleMemory[rankNo]));
    }
  }

  // deallocate mpi memory
  MPIUtility::handleReturnValue(MPI_Win_free(&mpiMemoryWindow), "MPI_Win_free");

  VLOG(1) << "after fence, nDofRequestedFromRanks: " << nDofRequestedFromRanks;

  // exchange which dofs are requested
  std::vector<MPI_Request> sendRequests;


  std::vector<std::vector<int>> sendBuffer(requestDofsFromRanks.size());

  int i = 0;
  for (typename std::map<int,DofsRequest>::iterator iter = requestDofsFromRanks.begin(); iter != requestDofsFromRanks.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nRequestedDofs = iter->second.dofNosLocal.size();

    sendBuffer[i].resize(nRequestedDofs);
    std::copy(iter->second.dofNosGlobalPetsc.begin(), iter->second.dofNosGlobalPetsc.end(), sendBuffer[i].begin());

    MPI_Request sendRequest;
    MPIUtility::handleReturnValue(MPI_Isend(sendBuffer[i].data(), nRequestedDofs, MPI_INT, foreignRankNo, 0,
                                            this->meshPartition_->mpiCommunicator(), &sendRequest), "MPI_Isend");

    VLOG(1) << "to rank " << foreignRankNo << " send " << nRequestedDofs << " requests: " << sendBuffer[i];

    sendRequests.push_back(sendRequest);
  }

  // receive which dofs are requested
  std::vector<std::vector<int>> requestedDofsGlobalPetsc(nDofRequestedFromRanks.size());
  std::vector<MPI_Request> receiveRequests;

  i = 0;
  for (typename std::vector<std::pair<int,int>>::iterator iter = nDofRequestedFromRanks.begin(); iter != nDofRequestedFromRanks.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nFromRank = iter->second;

    if (nFromRank != 0)
    {
      VLOG(1) << "i=" << i << ", from rank " << foreignRankNo << " receive " << nFromRank << " requests";

      requestedDofsGlobalPetsc[i].resize(nFromRank);
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(requestedDofsGlobalPetsc[i].data(), nFromRank, MPI_INT, foreignRankNo, 0,
                                              this->meshPartition_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);
    }
  }

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  sendRequests.clear();
  receiveRequests.clear();

  if (VLOG_IS_ON(1))
  {
    std::stringstream s;
    for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks)
    {
      s << "[rank " << requestDofsFromRank.first << ", dofNosGlobalPetsc: " << requestDofsFromRank.second.dofNosGlobalPetsc
        << ", dofNosLocal: " << requestDofsFromRank.second.dofNosLocal << "], ";
    }
    VLOG(1) << "requestDofsFromRanks: " << s.str();
    VLOG(1) << "requestedDofsGlobalPetsc: " << requestedDofsGlobalPetsc;
  }

  // send dofs and values in nonBcGlobal ordering
  i = 0;
  std::vector<std::vector<int>>    requestedDofsGlobalPetscSendBuffer(nDofRequestedFromRanks.size());
  std::vector<std::vector<double>> requestedDofsGlobalPetscSendBufferValues(nDofRequestedFromRanks.size());

  for (const std::pair<int,int> &nDofRequestedFromRank : nDofRequestedFromRanks)
  {
    int foreignRankNo = nDofRequestedFromRank.first;
    int nFromRank = nDofRequestedFromRank.second;

    if (nFromRank != 0)
    {
      requestedDofsGlobalPetscSendBuffer[i].resize(nFromRank*nComponentsDirichletBc);
      requestedDofsGlobalPetscSendBufferValues[i].resize(nFromRank*nComponentsDirichletBc, 0.0);

      for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
      {
        for (int requestedDofIndex = 0; requestedDofIndex < nFromRank; requestedDofIndex++)
        {
          global_no_t requestedDofNoGlobalPetsc = requestedDofsGlobalPetsc[i][requestedDofIndex];
          bool isLocal;
          dof_no_t requestedDofNoLocal = this->meshPartition_->getDofNoLocal(requestedDofNoGlobalPetsc, isLocal);

          requestedDofsGlobalPetscSendBuffer[i][componentNo*nFromRank + requestedDofIndex] = dofNoLocalToDofNoNonBcGlobal_[componentNo][requestedDofNoLocal];
          // dofNoLocalToDofNoNonBcGlobal_[componentNo][requestedDofNoLocal] is -1 if it is a dirichlet boundary condition dof

          requestedDofsGlobalPetscSendBufferValues[i][componentNo*nFromRank + requestedDofIndex] = this->boundaryConditionValues_[componentNo][requestedDofNoLocal];

          VLOG(1) << " send to rank " << foreignRankNo << ", component " << componentNo << " requested dofNoGlobalPetsc: " << requestedDofNoGlobalPetsc
            << ", requestedDofNoLocal: " << requestedDofNoLocal << " send dof " << dofNoLocalToDofNoNonBcGlobal_[componentNo][requestedDofNoLocal]
            << ", at sendBuffer index " << componentNo*nFromRank + requestedDofIndex;
        }
      }

      MPI_Request sendRequestDofs;
      MPIUtility::handleReturnValue(MPI_Isend(requestedDofsGlobalPetscSendBuffer[i].data(), nFromRank*nComponentsDirichletBc, MPI_INT, foreignRankNo, 0,
                                              this->meshPartition_->mpiCommunicator(), &sendRequestDofs), "MPI_Isend");
      sendRequests.push_back(sendRequestDofs);

      MPI_Request sendRequestValues;
      MPIUtility::handleReturnValue(MPI_Isend(requestedDofsGlobalPetscSendBufferValues[i].data(), nFromRank*nComponentsDirichletBc, MPI_DOUBLE, foreignRankNo, 1,
                                              this->meshPartition_->mpiCommunicator(), &sendRequestValues), "MPI_Isend");
      sendRequests.push_back(sendRequestValues);
    }
    i++;
  }

  // receive dofs and values in nonBCGlobal ordering
  std::vector<std::vector<int>> requestedDofsGlobalPetscReceiveBuffer(requestDofsFromRanks.size());
  std::vector<std::vector<double>> requestedDofsGlobalPetscReceiveBufferValues(requestDofsFromRanks.size());

  i = 0;
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks)
  {
    int foreignRankNo = requestDofsFromRank.first;
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();

    requestedDofsGlobalPetscReceiveBuffer[i].resize(nRequestedDofs*nComponentsDirichletBc);
    requestedDofsGlobalPetscReceiveBufferValues[i].resize(nRequestedDofs*nComponentsDirichletBc);

    MPI_Request receiveRequestDofs;
    MPIUtility::handleReturnValue(MPI_Irecv(requestedDofsGlobalPetscReceiveBuffer[i].data(), nRequestedDofs*nComponentsDirichletBc, MPI_INT, foreignRankNo, 0,
                                            this->meshPartition_->mpiCommunicator(), &receiveRequestDofs), "MPI_Irecv");
    receiveRequests.push_back(receiveRequestDofs);

    MPI_Request receiveRequestValues;
    MPIUtility::handleReturnValue(MPI_Irecv(requestedDofsGlobalPetscReceiveBufferValues[i].data(), nRequestedDofs*nComponentsDirichletBc, MPI_DOUBLE, foreignRankNo, 1,
                                            this->meshPartition_->mpiCommunicator(), &receiveRequestValues), "MPI_Irecv");
    receiveRequests.push_back(receiveRequestValues);
    i++;
  }

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  int dofNoNonBcLocal = nEntriesLocal_;

  // copy received dofs to new vector
  // and copy received values
  i = 0;
  nonBcGhostDofNosGlobal_.clear();
  VLOG(1) << "n recv buffer values: " << requestedDofsGlobalPetscReceiveBufferValues.size();
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks)
  {
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();

    //const std::vector<global_no_t> &ghostDofNosGlobalPetsc = requestDofsFromRank.second.dofNosGlobalPetsc;
    const std::vector<dof_no_t> &ghostDofNosLocal = requestDofsFromRank.second.dofNosLocal;

    for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
    {
      for (int j = 0; j < nRequestedDofs; j++)
      {
        global_no_t dofNoNonBcGlobal = requestedDofsGlobalPetscReceiveBuffer[i][componentNo*nRequestedDofs + j];
        double value = requestedDofsGlobalPetscReceiveBufferValues[i][componentNo*nRequestedDofs + j];
        dof_no_t dofNoLocal = ghostDofNosLocal[j];    // this value was not sent

        boundaryConditionValues_[componentNo][dofNoLocal] = value;

        // if the received dofNoNonBcGlobal does not contain the index, but the value -1, this means the dof has a Dirichlet BC
        if (dofNoNonBcGlobal == -1)
        {
          isPrescribed_[componentNo][dofNoLocal] = true;
          dofNoLocalToDofNoNonBcLocal_[componentNo][dofNoLocal] = -1;
        }
        else
        {
          nonBcGhostDofNosGlobal_.push_back(dofNoNonBcGlobal);
          dofNoLocalToDofNoNonBcLocal_[componentNo][dofNoLocal] = dofNoNonBcLocal;
          dofNoNonBcLocal++;
        }
        dofNoLocalToDofNoNonBcGlobal_[componentNo][dofNoLocal] = dofNoNonBcGlobal;

        VLOG(1) << "received from rank " << requestDofsFromRank.first << ", component " << componentNo << ", j: " << j << ", i: " << i
         << ", " << requestedDofsGlobalPetscReceiveBufferValues[i] << " entries, index=" << componentNo << "*" << nRequestedDofs << " + " << j << "=" << componentNo*nRequestedDofs + j
         << " dofNoLocal: " << dofNoLocal << ", value " << value << ", dofNoNonBcGlobal: " << dofNoNonBcGlobal << "(" << (dofNoNonBcGlobal==-1) << ")";

      }
    }
    i++;
  }
  VLOG(1) << "dofNoLocalToDofNoNonBcLocal_: " << dofNoLocalToDofNoNonBcLocal_;
  VLOG(1) << "dofNoLocalToDofNoNonBcGlobal_: " << dofNoLocalToDofNoNonBcGlobal_;
  VLOG(1) << "nonBcGhostDofNosGlobal_: " << nonBcGhostDofNosGlobal_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
createVector()
{
  nNonBcDofsGhosts_ = nonBcGhostDofNosGlobal_.size();

  VLOG(1) << "nEntriesLocal_: " << nEntriesLocal_ << ", nEntriesGlobal_: " << nEntriesGlobal_ << ", nNonBcDofsGhosts_: " << nNonBcDofsGhosts_;
  VLOG(1) << "nonBcGhostDofNosGlobal_: " << nonBcGhostDofNosGlobal_;
  VLOG(1) << "boundaryConditionValues_: " << boundaryConditionValues_ << ", isPrescribed_: " << isPrescribed_;

  // initialize PETSc vector object
  PetscErrorCode ierr;
  ierr = VecCreateGhost(this->meshPartition_->mpiCommunicator(), nEntriesLocal_,
                        nEntriesGlobal_, nNonBcDofsGhosts_, nonBcGhostDofNosGlobal_.data(), &vectorCombinedWithoutDirichletDofsGlobal_); CHKERRV(ierr);

  ierr = PetscObjectSetName((PetscObject) vectorCombinedWithoutDirichletDofsGlobal_, this->name_.c_str()); CHKERRV(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(vectorCombinedWithoutDirichletDofsGlobal_); CHKERRV(ierr);
  ierr = VecGhostGetLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

  ierr = VecSetOption(vectorCombinedWithoutDirichletDofsLocal_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE); CHKERRV(ierr);
  ierr = VecSetOption(vectorCombinedWithoutDirichletDofsGlobal_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE); CHKERRV(ierr);

  // createVector acts like startGhostManipulation as it also gets the local vector (VecGhostGetLocalForm) to work on.
  this->currentRepresentation_ = Partition::values_representation_t::representationCombinedLocal;
}

//! get number of local entries
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
int PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nEntriesLocal() const
{
  return nEntriesLocal_;
}

//! get number of global entries
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
global_no_t PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
nEntriesGlobal() const
{
  return nEntriesGlobal_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
Vec &PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
valuesGlobal()
{
  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedGlobal)
  {
    VLOG(1) << "valuesGlobal called in not combined-global vector representation ("
      << Partition::valuesRepresentationString[this->currentRepresentation_]
      <<"), now set to combined-global (without considering ghost dofs, call finishGhostManipulation if ghosts are needed!)";
    setRepresentationGlobal();
  }

  VLOG(2) << "\"" << this->name_ << "\" valuesGlobal()";

  LOG(DEBUG) << "valuesGlobal, return vectorCombinedWithoutDirichletDofsGlobal_";

  return vectorCombinedWithoutDirichletDofsGlobal_;
}

//! Communicates the ghost values from the global vectors to the local vector and sets the representation to local.
//! The representation has to be global, afterwards it is set to local.
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
startGhostManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" startGhostManipulation";

  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedGlobal)
  {
    LOG(FATAL) << "\"" << this->name_ << "\", startGhostManipulation called when representation is not combined-global (but "
      << Partition::valuesRepresentationString[this->currentRepresentation_]
      << "), this overwrites the previous values and fetches the last from the global vectors!" << std::endl
      << "Call setRepresentationGlobal() before startGhostManipulation() or check if startGhostManipulation() "
      << "is even necessary (because the representation is already local).";
  }

  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;

  ierr = VecGhostUpdateBegin(vectorCombinedWithoutDirichletDofsGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
  ierr = VecGhostUpdateEnd(vectorCombinedWithoutDirichletDofsGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
  ierr = VecGhostGetLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

  this->currentRepresentation_ = Partition::values_representation_t::representationCombinedLocal;
}

//! Communicates the ghost values from the local vectors back to the global vector and sets the representation to global.
//! The representation has to be local, afterwards it is set to global.
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
finishGhostManipulation()
{
  VLOG(2) << "\"" << this->name_ << "\" finishGhostManipulation";

  if (this->currentRepresentation_ != Partition::values_representation_t::representationCombinedLocal)
  {
    LOG(ERROR) << "\"" << this->name_ << "\", finishGhostManipulation called when representation is not combined-local (it is "
      << Partition::valuesRepresentationString[this->currentRepresentation_]
      << "), (probably no previous startGhostManipulation)";
  }

  // Copy the local values vectors into the global vector. ADD_VALUES means that ghost values are reduced (summed up)
  PetscErrorCode ierr;
  ierr = VecGhostRestoreLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);
  ierr = VecGhostUpdateBegin(vectorCombinedWithoutDirichletDofsGlobal_, ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);
  ierr = VecGhostUpdateEnd(vectorCombinedWithoutDirichletDofsGlobal_, ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);

  this->currentRepresentation_ = Partition::values_representation_t::representationCombinedGlobal;
}

// set the internal representation to be global, i.e. using the global vectors
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setRepresentationGlobal()
{
  VLOG(2) << "\"" << this->name_ << "\" PartitionedPetscVecWithDirichletBc::setRepresentationGlobal, previous representation: "
    << Partition::valuesRepresentationString[this->currentRepresentation_];

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    return;
  }

    // set the internal representation to be global, i.e. using the global vectors
  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    // restore local vector back into global vector, does not consider ghost dofs
    // use finishGhostManipulation to handle ghost dofs, that requires communication

    PetscErrorCode ierr;
    ierr = VecGhostRestoreLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

    this->currentRepresentation_ = Partition::values_representation_t::representationCombinedGlobal;
  }
  else
  {
    LOG(FATAL) << "Cannot set vector representation from \"" << Partition::valuesRepresentationString[this->currentRepresentation_]
      << "\" to \"combined-global\".";
  }
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setRepresentationLocal()
{
  VLOG(2) << "\"" << this->name_ << "\" setRepresentationLocal, previous representation: "
    << Partition::valuesRepresentationString[this->currentRepresentation_];

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    // retrieve the local vector from the global vector, does not consider ghost dofs
    // use startGhostManipulation to also set the correct values for the ghost dofs, that requires communication

    PetscErrorCode ierr;
    ierr = VecGhostGetLocalForm(vectorCombinedWithoutDirichletDofsGlobal_, &vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);

    this->currentRepresentation_ = Partition::values_representation_t::representationCombinedLocal;
  }
  else
  {
    LOG(FATAL) << "Cannot set vector representation from \"" << Partition::valuesRepresentationString[this->currentRepresentation_]
      << "\" to \"combined-local\".";
  }
}


//! zero all values in the local ghost buffer. Needed if between startGhostManipulation() and finishGhostManipulation() only some ghost will be reassigned. To prevent that the "old" ghost values that were present in the local ghost values buffer get again added to the real values which actually did not change.
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
zeroGhostBuffer()
{
  VLOG(2) << "\"" << this->name_ << "\" zeroGhostBuffer";

  // set local ghost values to 0
  std::vector<int> indices(nNonBcDofsGhosts_);
  std::iota(indices.begin(), indices.end(), nEntriesLocal_);
  std::vector<double> values(nNonBcDofsGhosts_, 0.0);

  PetscErrorCode ierr;
  ierr = VecSetValues(vectorCombinedWithoutDirichletDofsLocal_, nNonBcDofsGhosts_, indices.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
}

//! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
{
  VLOG(3) << "\"" << this->name_ << "\" setValues, representation: " << Partition::valuesRepresentationString[this->currentRepresentation_];

  assert(componentNo >= 0 && componentNo < nComponents);

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    // determine new indices
    std::vector<int> indices;
    std::vector<double> values;

    indices.reserve(ni);
    values.reserve(ni);

  PetscErrorCode ierr;
#ifndef NDEBUG
    int ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vectorCombinedWithoutDirichletDofsGlobal_, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
#endif

    for (int i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithoutGhosts());

      if (!isPrescribed_[componentNo][ix[i]])
      {
        int nonBcIndexGlobal = nonBCDofNoGlobal(componentNo, ix[i]);
#ifndef NDEBUG
        if (!(nonBcIndexGlobal >= ownershipBegin && nonBcIndexGlobal < ownershipEnd))
        {
          LOG(ERROR) << "setValues \"" << this->name_ << "\", i=" << i << ", ix[i]=" << ix[i] << ", nonBcIndexGlobal=" << nonBcIndexGlobal
            << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
        }
        assert(nonBcIndexGlobal >= ownershipBegin && nonBcIndexGlobal < ownershipEnd);
#endif

        indices.push_back(nonBcIndexGlobal);
        values.push_back(y[i]);
      }
    }

    VLOG(1) << "non-bc indices: " << indices << ", values: " << values;
    // this wraps the standard PETSc VecSetValues on the local vector
    ierr = VecSetValues(vectorCombinedWithoutDirichletDofsGlobal_, indices.size(), indices.data(), values.data(), iora); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
  else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    // determine new indices
    std::vector<int> indices;
    std::vector<double> values;

    indices.reserve(ni);
    values.reserve(ni);

    for (int i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithGhosts());

      if (!isPrescribed_[componentNo][ix[i]])
      {
        indices.push_back(nonBCDofNoLocal(componentNo, ix[i]));
        values.push_back(y[i]);
      }
    }

    VLOG(1) << "non-bc indices: " << indices << ", values: " << values;
    // this wraps the standard PETSc VecSetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValues(vectorCombinedWithoutDirichletDofsLocal_, indices.size(), indices.data(), values.data(), iora); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
}

//! wrapper to the PETSc VecSetValue, acting only on the local data
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode)
{
  VLOG(3) << "\"" << this->name_ << "\" setValue row=" << row << ", representation: " << Partition::valuesRepresentationString[this->currentRepresentation_];

  assert(componentNo >= 0 && componentNo < nComponents);

  // replace dirichlet BC values with the prescribed values
  if (isPrescribed_[componentNo][row])
  {
    VLOG(1) << "row " << row << ", value is prescribed, do not change";
    return;
  }

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    assert(row < this->meshPartition_->nDofsLocalWithoutGhosts());

  PetscErrorCode ierr;
#ifndef NDEBUG
    int ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vectorCombinedWithoutDirichletDofsGlobal_, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
#endif

    // determine new index
    row = nonBCDofNoGlobal(componentNo, row);

#ifndef NDEBUG
    if (!(row >= ownershipBegin && row < ownershipEnd))
    {
      LOG(ERROR) << "setValue \"" << this->name_ << "\", row=" << row
        << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
    }
    assert(row >= ownershipBegin && row < ownershipEnd);
#endif

    VLOG(1) << "set value " << value << " at non-bc global " << row;

    // this wraps the standard PETSc VecSetValues on the local vector
    ierr = VecSetValue(vectorCombinedWithoutDirichletDofsGlobal_, row, value, mode); CHKERRV(ierr);
  }
  else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    assert(row < this->meshPartition_->nDofsLocalWithGhosts());

    // determine new index
    row = nonBCDofNoLocal(componentNo, row);

    VLOG(1) << "set value " << value << " at non-bc local " << row;

    // this wraps the standard PETSc VecSetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValue(vectorCombinedWithoutDirichletDofsLocal_, row, value, mode); CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but the indices must be provided in the local ordering.
    // For this to work one has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
}

//! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const
{
  VLOG(3) << "\"" << this->name_ << "\" getValues, representation: " << Partition::valuesRepresentationString[this->currentRepresentation_];

  assert(componentNo >= 0 && componentNo < nComponents);

  if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {

    PetscErrorCode ierr;
#ifndef NDEBUG
    int ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vectorCombinedWithoutDirichletDofsGlobal_, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);

    // check that indices are in range
    for (int i = 0; i < ni; i++)
    {
      if (ix[i] >= this->meshPartition_->nDofsLocalWithoutGhosts())
      {
        LOG(ERROR) << "ni = " << ni << ", ix[" << i << "]=" << ix[i] << " >= " << this->meshPartition_->nDofsLocalWithoutGhosts();
      }
    }
#endif

    // determine new global indices
    std::vector<int> indices(ni);
    for (int i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithoutGhosts());
      indices[i] = nonBCDofNoGlobal(componentNo, ix[i]);

      // check if indices are in range that is owned by own rank
#ifndef NDEBUG
      if (!(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd)))
      {
        LOG(ERROR) << "getValues \"" << this->name_ << "\", i=" << i << ", ix[i]=" << ix[i] << ", indices[i]=" << indices[i]
          << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
      }
      assert(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd));
#endif
    }

    VLOG(1) << "non-bc global indices: " << indices;

    // this wraps the standard PETSc VecGetValues on the local vector
    ierr = VecGetValues(vectorCombinedWithoutDirichletDofsGlobal_, ni, indices.data(), y); CHKERRV(ierr);

    VLOG(1) << "got values: ";
    for (int i = 0; i < ni; i++)
      VLOG(1) << y[i];

    // replace dirichlet BC values with the prescribed values
    for (int i = 0; i < ni; i++)
    {
      if (isPrescribed_[componentNo][ix[i]] && componentNo < nComponentsDirichletBc)
      {
        y[i] = boundaryConditionValues_[componentNo][ix[i]];
      }
    }

    VLOG(1) << "modified values: ";
    for (int i = 0; i < ni; i++)
      VLOG(1) << y[i];
  }
  else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    // determine new indices
    std::vector<int> indices(ni);
    for (int i = 0; i < ni; i++)
    {
      assert(ix[i] < this->meshPartition_->nDofsLocalWithGhosts());
      indices[i] = nonBCDofNoLocal(componentNo, ix[i]);
    }

    VLOG(1) << "non-bc local indices: " << indices;

    // this wraps the standard PETSc VecGetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecGetValues(vectorCombinedWithoutDirichletDofsLocal_, ni, indices.data(), y); CHKERRV(ierr);

    VLOG(1) << "got values: ";
    for (int i = 0; i < ni; i++)
      VLOG(1) << y[i];

    // replace dirichlet BC values with the prescribed values
    for (int i = 0; i < ni; i++)
    {
      if (isPrescribed_[componentNo][ix[i]] && componentNo < nComponentsDirichletBc)
      {
        y[i] = boundaryConditionValues_[componentNo][ix[i]];
      }
    }

    VLOG(1) << "modified values: ";
    for (int i = 0; i < ni; i++)
      VLOG(1) << y[i];
  }
}

//! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
getValuesGlobal(Vec vector, int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const
{
  LOG(FATAL) << "getValuesGlobal should not be called.";

  PetscErrorCode ierr;
  std::string name;
  char *cName;
  ierr = PetscObjectGetName((PetscObject)vector, (const char **)&cName); CHKERRV(ierr);
  name = cName;

  VLOG(3) << " getValuesGlobal(\"" << name << "\")";

#ifndef NDEBUG
    int ownershipBegin, ownershipEnd;
    ierr = VecGetOwnershipRange(vector, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
#endif

  // determine new global indices
  std::vector<int> indices(ni);
  for (int i = 0; i < ni; i++)
  {
    assert(ix[i] < this->meshPartition_->nDofsLocalWithoutGhosts());
    indices[i] = nonBCDofNoGlobal(componentNo, ix[i]);

#ifndef NDEBUG
    if (!(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd)))
    {
      LOG(ERROR) << "getValuesGlobal \"" << this->name_ << "\", i=" << i << ", ix[i]=" << ix[i] << ", indices[i]=" << indices[i]
        << ", ownership: [" << ownershipBegin << "," << ownershipEnd << "]";
    }
    assert(indices[i] == -1 || (indices[i] >= ownershipBegin && indices[i] < ownershipEnd));
#endif
  }

  VLOG(1) << "non-bc global indices: " << indices;

  // this wraps the standard PETSc VecGetValues on the local vector
  ierr = VecGetValues(vector, ni, indices.data(), y); CHKERRV(ierr);

  VLOG(1) << "got values: ";
  for (int i = 0; i < ni; i++)
    VLOG(1) << y[i];

  // replace dirichlet BC values with the prescribed values
  for (int i = 0; i < ni; i++)
  {
    if (isPrescribed_[componentNo][ix[i]])
    {
      if (componentNo < nComponentsDirichletBc)
      {
        y[i] = boundaryConditionValues_[componentNo][ix[i]];
      }
    }
  }

  VLOG(1) << "modified values: ";
  for (int i = 0; i < ni; i++)
    VLOG(1) << y[i];
}

//! get a single value
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
double PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
getValue(int componentNo, PetscInt row) const
{
  double result;
  this->getValues(componentNo, 1, &row, &result);
  return result;
}

//! set all entries to zero, wraps VecZeroEntries
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
zeroEntries()
{
  VLOG(3) << "\"" << this->name_ << "\" zeroEntries, representation: " << Partition::valuesRepresentationString[this->currentRepresentation_];

  PetscErrorCode ierr;
  //if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedLocal)
  {
    ierr = VecZeroEntries(vectorCombinedWithoutDirichletDofsLocal_); CHKERRV(ierr);
  }
  //else if (this->currentRepresentation_ == Partition::values_representation_t::representationCombinedGlobal)
  {
    ierr = VecZeroEntries(vectorCombinedWithoutDirichletDofsGlobal_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
updateDirichletBoundaryConditions(const std::vector<std::pair<global_no_t,std::array<double,nComponentsDirichletBc>>> &newValues,
                                  bool inputMeshIsGlobal)
{
  for (int itemIndex = 0; itemIndex < newValues.size(); itemIndex++)
  {
    dof_no_t dofNoLocal = newValues[itemIndex].first;

    if (inputMeshIsGlobal)
    {
      global_no_t dofNoGlobal = newValues[itemIndex].first;

      bool isLocal = false;
      dofNoLocal = this->meshPartition_->getDofNoLocal(dofNoGlobal, isLocal);

      if (!isLocal)
        continue;
    }

    std::array<double,nComponentsDirichletBc> values = newValues[itemIndex].second;

    for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
    {
      boundaryConditionValues_[componentNo][dofNoLocal] = values[componentNo];
    }
  }
}


//! output the vector to stream, for debugging
template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
output(std::ostream &stream) const
{
  stream << "\"" << this->name_ << "\" " << this->nEntriesGlobal_ << " entries global, " << nEntriesLocal_ << " entries local." << std::endl
    << "nNonBcDofsWithoutGhosts: " << nNonBcDofsWithoutGhosts_ << ", nNonBcDofsGhosts: " << nNonBcDofsGhosts_ << ", nonBcGhostDofNosGlobal_: " << nonBcGhostDofNosGlobal_ << std::endl
    << "dofNoLocalToDofNoNonBcGlobal_: " << dofNoLocalToDofNoNonBcGlobal_ << std::endl
    << "dofNoLocalToDofNoNonBcLocal_:  " << dofNoLocalToDofNoNonBcLocal_ << std::endl
    << "boundaryConditionValues_: " << boundaryConditionValues_ << ", isPrescribed_: " << isPrescribed_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscVecWithDirichletBc<FunctionSpaceType,nComponents> &vector)
{
  vector.output(stream);
  return stream;
}
