#include "partition/partitioned_petsc_vec/01_partitioned_petsc_vec_with_dirichlet_bc.h"

#include "utility/mpi_utility.h"
#include "spatial_discretization/dirichlet_boundary_conditions/00_dirichlet_boundary_conditions_base.h"

#define USE_MPI_RMA       // whether to use MPI rma, if not, use MPI_Alltoall

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
initialize(int offsetInGlobalNumberingPerRank)
{
  LOG(DEBUG) << "\"" << this->name_ << "\" PartitionedPetscVecWithDirichletBc createVector with " << nComponents << "," << nComponentsDirichletBc 
    << " components, size local: " << this->meshPartition_->nNodesLocalWithoutGhosts()
    << ", global: " << this->meshPartition_->nNodesGlobal() << ", offsetInGlobalNumberingPerRank: " << offsetInGlobalNumberingPerRank
    << ", ghost dof nos global/petsc: " << this->meshPartition_->ghostDofNosGlobalPetsc();

  VLOG(1) << "meshPartition: " << *this->meshPartition_;

  if (!nDofRequestedFromRanks_.empty())
  {
    LOG(FATAL) << "is already initialized";
  }

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

  // add space for pressure components, this is only !=0 for PartitionedPetscVecForHyperelasticity for incompressible formulation
  nDofsLocal_ = nEntriesLocal_;
  nEntriesLocal_ += offsetInGlobalNumberingPerRank;

  LOG(DEBUG) << "nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts << ", nDofsLocal_: " << nDofsLocal_ << ", nEntriesLocal_: " << nEntriesLocal_;

  // determine global number of non-BC dofs over all components
  nEntriesGlobal_ = 0;
  MPIUtility::handleReturnValue(MPI_Allreduce(&nEntriesLocal_, &nEntriesGlobal_, 1, MPI_INT, MPI_SUM, this->meshPartition_->mpiCommunicator()), "MPI_Allreduce");

  // get number of non-bc dofs on previous ranks
  nonBcDofNoGlobalBegin_ = 0;
  MPIUtility::handleReturnValue(MPI_Exscan(&nEntriesLocal_, &nonBcDofNoGlobalBegin_, 1, MPI_INT, MPI_SUM, this->meshPartition_->mpiCommunicator()), "MPI_Exscan");

  // setup non-bc numberings (local and global) for non-ghost dofs, also store boundaryConditionValues_ and isPrescribed_ for non-ghost dofs
  global_no_t dofNoNonBcGlobal = nonBcDofNoGlobalBegin_;

  VLOG(1) << "nEntriesLocal_: " << nEntriesLocal_ << ", nEntriesGlobal_: " << nEntriesGlobal_ << ", nonBcDofNoGlobalBegin_: " << nonBcDofNoGlobalBegin_;
  
  // loop over components
  for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
  {
    std::vector<dof_no_t>::const_iterator boundaryConditionDofLocalNosIter = dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].dofNosLocal.begin();
    std::vector<double>::const_iterator boundaryConditionValuesIter = dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].values.begin();

    LOG(DEBUG) << "number of boundary conditions for component " << componentNo << ": " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].values.size();

    dofNoLocalToDofNoNonBcGlobal_[componentNo].resize(nDofsLocalWithGhosts);
    dofNoLocalToDofNoNonBcLocal_[componentNo].resize(nDofsLocalWithGhosts);
    boundaryConditionValues_[componentNo].resize(nDofsLocalWithGhosts, -2);  // these initial values should all be overwritten later, by -1 for not prescribed values or by the actual value
    isPrescribed_[componentNo].resize(nDofsLocalWithGhosts, false);

    // loop over local dofs
    for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
    {
      VLOG(2) << "componentNo " << componentNo << " dofNoLocal " << dofNoLocal;

      // if dof has bc
      if (boundaryConditionDofLocalNosIter != dirichletBoundaryConditions_->boundaryConditionsByComponent()[componentNo].dofNosLocal.end())
      {
        VLOG(2) << " next bc: " << *boundaryConditionDofLocalNosIter;

        if (*boundaryConditionDofLocalNosIter == dofNoLocal)
        {
          VLOG(2) << " -> is bc";

          dofNoLocalToDofNoNonBcGlobal_[componentNo][dofNoLocal] = -1;
          dofNoLocalToDofNoNonBcLocal_[componentNo][dofNoLocal] = -1;
          boundaryConditionValues_[componentNo][dofNoLocal] = *boundaryConditionValuesIter;
          isPrescribed_[componentNo][dofNoLocal] = true;

          boundaryConditionDofLocalNosIter++;
          boundaryConditionValuesIter++;
          continue;
        }
      }
      VLOG(2) << " is no bc, assign new no. " << dofNoNonBcGlobal;

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
/**
  * struct DofsRequest
  *  {
  *    std::vector<global_no_t> dofNosGlobalPetsc;
  *    std::vector<dof_no_t> dofNosLocal;
  *  };
  *
  *  std::map<int, DofsRequest> requestDofsFromRanks_;
  */

  LOG(DEBUG) << "loop over ghost nodes";
  for (dof_no_t dofNoLocal = nDofsLocalWithoutGhosts; dofNoLocal < this->meshPartition_->nDofsLocalWithGhosts(); dofNoLocal++)
  {
    global_no_t dofNoGlobalPetsc = this->meshPartition_->getDofNoGlobalPetsc(dofNoLocal);

    // if the current dof is a ghost dof on the own rank (should be)
    node_no_t nodeNoLocal = dofNoLocal / FunctionSpaceType::nDofsPerNode();

    VLOG(1) << "dofNoLocal: " << dofNoLocal << "/[" << nDofsLocalWithoutGhosts << "," << this->meshPartition_->nDofsLocalWithGhosts()
      << "] nodeNoLocal: " << nodeNoLocal << ", dofNoGlobalPetsc: " << dofNoGlobalPetsc;

    int neighbourRankNo = 0;
    if (!this->meshPartition_->isNonGhost(nodeNoLocal, neighbourRankNo))
    {
      VLOG(1) << "request this dof from rank " << neighbourRankNo << " (this was determined by this->meshPartition_->isNonGhost)";
      requestDofsFromRanks_[neighbourRankNo].dofNosLocal.push_back(dofNoLocal);
      requestDofsFromRanks_[neighbourRankNo].dofNosGlobalPetsc.push_back(dofNoGlobalPetsc);
    }
    else
    {
      LOG(DEBUG) << "isNonGhost returned true, neighbourRankNo: " << neighbourRankNo;
      LOG(FATAL) << "ghost dof not recognized as ghost dof (isNonGhost is errorneous)";
    }
  }

  // exchange, how many values should be sent to which rank
  VLOG(1) << "rankSubset " << *this->meshPartition_->rankSubset() << ", create new window";

  int nRanks = this->meshPartition_->nRanks();
  int ownRankNo = this->meshPartition_->ownRankNo();

#ifdef USE_MPI_RMA   // implementation using RMA, not very stable with OpenMPI
  // create remote accessible memory
  MPI_Win mpiMemoryWindow;
  int displacementUnit = sizeof(int);
  int nBytes = nRanks * sizeof(int);

#ifdef USE_MPI_ALLOC
  LOG(DEBUG) << "Using MPI_Win_allocate";

  // let MPI allocate the memory
  int *remoteAccessibleMemory = nullptr;
  MPIUtility::handleReturnValue(MPI_Win_allocate(nBytes, displacementUnit, MPI_INFO_NULL, this->meshPartition_->mpiCommunicator(), (void *)&remoteAccessibleMemory, &mpiMemoryWindow), "MPI_Win_allocate");

  // clear buffer
  memset(remoteAccessibleMemory, 0, nBytes);
#else
  LOG(DEBUG) << "Using MPI_Win_create";

  // allocate the memory ourselves
  std::vector<int> remoteAccessibleMemory(nRanks, 0);
  MPIUtility::handleReturnValue(MPI_Win_create((void *)remoteAccessibleMemory.data(), nBytes, displacementUnit, MPI_INFO_NULL, this->meshPartition_->mpiCommunicator(), &mpiMemoryWindow), "MPI_Win_create");
#endif

#endif

  std::vector<int> localMemory(nRanks);

  // fill send buffer in localMemory
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
  {
    int foreignRankNo = requestDofsFromRank.first;
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();
    localMemory[foreignRankNo] = nRequestedDofs;
  }


#ifdef USE_MPI_RMA   // implementation using RMA, not very stable with OpenMPI
  // put number of requested ghost dofs to the corresponding processes
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
  {
    int foreignRankNo = requestDofsFromRank.first;

    VLOG(1) << "put value " << localMemory[foreignRankNo] << " to rank " << foreignRankNo << ", offset " << ownRankNo;

    // start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
    MPIUtility::handleReturnValue(MPI_Win_lock(MPI_LOCK_SHARED, foreignRankNo, 0, mpiMemoryWindow), "MPI_Win_lock");
    //MPIUtility::handleReturnValue(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, foreignRankNo, 0, mpiMemoryWindow), "MPI_Win_lock");

    MPIUtility::handleReturnValue(MPI_Put(localMemory.data()+foreignRankNo, 1, MPI_INT, foreignRankNo, ownRankNo, 1, MPI_INT, mpiMemoryWindow), "MPI_Put");

    MPIUtility::handleReturnValue(MPI_Win_unlock(foreignRankNo, mpiMemoryWindow), "MPI_Win_unlock");
  }

  MPIUtility::handleReturnValue(MPI_Win_fence(MPI_MODE_NOSUCCEED, mpiMemoryWindow), "MPI_Win_fence");
#else
  LOG(DEBUG) << "Using MPI_Alltoall with rankSubset " << *this->meshPartition_->rankSubset() << ", local data: " << localMemory;

  // allocate the receiving memory
  std::vector<int> remoteAccessibleMemory(nRanks, 0);

  // alternative implementation using Alltoall
  MPIUtility::handleReturnValue(MPI_Alltoall(localMemory.data(), 1, MPI_INT, remoteAccessibleMemory.data(), 1, MPI_INT,
                                             this->meshPartition_->rankSubset()->mpiCommunicator()), "MPI_Alltoall");

  LOG(DEBUG) << "MPI_Alltoall is done, result: " << remoteAccessibleMemory;
#endif

  //std::vector<std::pair<int,int>> nDofRequestedFromRanks_;   /// (foreignRank,nDofs), number of dofs requested by and to be send to foreignRank
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    VLOG(1) << " rank " << rankNo << " nRequestedDofs: " << remoteAccessibleMemory[rankNo];
    if (remoteAccessibleMemory[rankNo] > 0)
    {
      nDofRequestedFromRanks_.push_back(std::pair<int,int>(rankNo,remoteAccessibleMemory[rankNo]));
    }
  }

#ifdef USE_MPI_RMA
#ifdef USE_MPI_ALLOC
  // deallocate mpi memory
  MPIUtility::handleReturnValue(MPI_Win_free(&mpiMemoryWindow), "MPI_Win_free");
#endif
#endif

  VLOG(1) << "after fence, nDofRequestedFromRanks_: " << nDofRequestedFromRanks_;

  // exchange which dofs are requested
  std::vector<MPI_Request> sendRequests;


  std::vector<std::vector<int>> sendBuffer(requestDofsFromRanks_.size());

  int i = 0;
  for (typename std::map<int,DofsRequest>::iterator iter = requestDofsFromRanks_.begin(); iter != requestDofsFromRanks_.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nRequestedDofs = iter->second.dofNosLocal.size();

    sendBuffer[i].resize(nRequestedDofs);
    std::copy(iter->second.dofNosGlobalPetsc.begin(), iter->second.dofNosGlobalPetsc.end(), sendBuffer[i].begin());
    assert(iter->second.dofNosGlobalPetsc.size() == nRequestedDofs);

    int tag = foreignRankNo*10000+nRequestedDofs;
    MPI_Request sendRequest;
    MPIUtility::handleReturnValue(MPI_Isend(sendBuffer[i].data(), nRequestedDofs, MPI_INT, foreignRankNo, tag,
                                            this->meshPartition_->mpiCommunicator(), &sendRequest), "MPI_Isend");

    VLOG(1) << "to rank " << foreignRankNo << " send " << nRequestedDofs << " requests: " << sendBuffer[i] << ", tag=" << tag;

    sendRequests.push_back(sendRequest);
  }

  // receive which dofs are requested
  requestedDofsGlobalPetsc_.resize(nDofRequestedFromRanks_.size());
  std::vector<MPI_Request> receiveRequests;

  i = 0;
  for (typename std::vector<std::pair<int,int>>::iterator iter = nDofRequestedFromRanks_.begin(); iter != nDofRequestedFromRanks_.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nFromRank = iter->second;

    if (nFromRank != 0)
    {
      int tag = ownRankNo*10000+nFromRank;
      VLOG(1) << "i=" << i << ", from rank " << foreignRankNo << " receive " << nFromRank << " requests, tag=" << tag;

      requestedDofsGlobalPetsc_[i].resize(nFromRank, -1);
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(requestedDofsGlobalPetsc_[i].data(), nFromRank, MPI_INT, foreignRankNo, tag,
                                              this->meshPartition_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);
    }
  }

  VLOG(1) << "MPI_Waitall (1), " << sendRequests.size() << " sendRequests, " << *this->meshPartition_->rankSubset();

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  VLOG(1) << "MPI_Waitall (1), " << receiveRequests.size() << " receiveRequests, " << *this->meshPartition_->rankSubset();

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  sendRequests.clear();
  receiveRequests.clear();

  if (VLOG_IS_ON(1))
  {
    std::stringstream s;
    for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
    {
      s << "[rank " << requestDofsFromRank.first << ", dofNosGlobalPetsc: " << requestDofsFromRank.second.dofNosGlobalPetsc
        << ", dofNosLocal: " << requestDofsFromRank.second.dofNosLocal << "], ";
    }
    VLOG(1) << "data to receive, as requested: requestDofsFromRanks_: " << s.str();
    VLOG(1) << "data to send to fulfill requests: requestedDofsGlobalPetsc_: " << requestedDofsGlobalPetsc_;
  }

  // send dofs and values in nonBcGlobal ordering
  i = 0;
  std::vector<std::vector<int>>    requestedDofsGlobalPetscSendBuffer(nDofRequestedFromRanks_.size());
  std::vector<std::vector<double>> requestedDofsGlobalPetscSendBufferValues(nDofRequestedFromRanks_.size());

  for (const std::pair<int,int> &nDofRequestedFromRank : nDofRequestedFromRanks_)
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
          global_no_t requestedDofNoGlobalPetsc = requestedDofsGlobalPetsc_[i][requestedDofIndex];
          bool isLocal;
          dof_no_t requestedDofNoLocal = this->meshPartition_->getDofNoLocal(requestedDofNoGlobalPetsc, isLocal);

          if (!isLocal)
          {
            LOG(FATAL) << "Error in partitioning, global dof no " << requestedDofNoGlobalPetsc << " is not on rank " << ownRankNo << ", but rank " << foreignRankNo << " thought it would be.";
          }

          requestedDofsGlobalPetscSendBuffer[i][componentNo*nFromRank + requestedDofIndex] = dofNoLocalToDofNoNonBcGlobal_[componentNo][requestedDofNoLocal];
          // dofNoLocalToDofNoNonBcGlobal_[componentNo][requestedDofNoLocal] is -1 if it is a dirichlet boundary condition dof

          requestedDofsGlobalPetscSendBufferValues[i][componentNo*nFromRank + requestedDofIndex] = this->boundaryConditionValues_[componentNo][requestedDofNoLocal];

          VLOG(1) << " send to rank " << foreignRankNo << ", component " << componentNo << " requested dofNoGlobalPetsc: " << requestedDofNoGlobalPetsc
            << ", requestedDofNoLocal: " << requestedDofNoLocal << " (isLocal: " << isLocal << " should be true), send dof " << dofNoLocalToDofNoNonBcGlobal_[componentNo][requestedDofNoLocal]
            << " and bc value " << this->boundaryConditionValues_[componentNo][requestedDofNoLocal] << " at sendBuffer index " << componentNo*nFromRank + requestedDofIndex;
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
  std::vector<std::vector<int>> requestedDofsGlobalPetscReceiveBuffer(requestDofsFromRanks_.size());
  std::vector<std::vector<double>> requestedDofsGlobalPetscReceiveBufferValues(requestDofsFromRanks_.size());

  i = 0;
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
  {
    int foreignRankNo = requestDofsFromRank.first;
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();

    requestedDofsGlobalPetscReceiveBuffer[i].resize(nRequestedDofs*nComponentsDirichletBc);
    requestedDofsGlobalPetscReceiveBufferValues[i].resize(nRequestedDofs*nComponentsDirichletBc);

    VLOG(1) << " receive from rank " << foreignRankNo << ", " << nRequestedDofs << " dofs and values";

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

  VLOG(1) << "MPI_Waitall (2), " << sendRequests.size() << " sendRequests, " << *this->meshPartition_->rankSubset();

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  VLOG(1) << "MPI_Waitall (2), " << receiveRequests.size() << " receiveRequests, " << *this->meshPartition_->rankSubset();

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  int dofNoNonBcLocal = nEntriesLocal_;

  // copy received dofs to new vector
  // and copy received values
  i = 0;
  nonBcGhostDofNosGlobal_.clear();
  VLOG(1) << "n recv buffer values: " << requestedDofsGlobalPetscReceiveBufferValues.size();
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
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
         << ", " << requestedDofsGlobalPetscReceiveBufferValues[i].size() << " entries, index=" << componentNo << "*" << nRequestedDofs << " + " << j << "=" << componentNo*nRequestedDofs + j
         << " dofNoLocal: " << dofNoLocal << ", value " << value << ", dofNoNonBcGlobal: " << dofNoNonBcGlobal << "(" << (dofNoNonBcGlobal==-1) << ")";

      }
    }
    i++;
  }

  nNonBcDofsGhosts_ = nonBcGhostDofNosGlobal_.size();

  VLOG(1) << "dofNoLocalToDofNoNonBcLocal_: " << dofNoLocalToDofNoNonBcLocal_;
  VLOG(1) << "dofNoLocalToDofNoNonBcGlobal_: " << dofNoLocalToDofNoNonBcGlobal_;
  VLOG(1) << "nonBcGhostDofNosGlobal_: " << nonBcGhostDofNosGlobal_;
}

template<typename FunctionSpaceType, int nComponents, int nComponentsDirichletBc>
void PartitionedPetscVecWithDirichletBc<FunctionSpaceType, nComponents, nComponentsDirichletBc>::
communicateBoundaryConditionGhostValues()
{
  std::vector<MPI_Request> receiveRequests;
  std::vector<MPI_Request> sendRequests;


  // send dofs and values in nonBcGlobal ordering
  int i = 0;
  std::vector<std::vector<double>> requestedDofsGlobalPetscSendBufferValues(nDofRequestedFromRanks_.size());

  for (const std::pair<int,int> &nDofRequestedFromRank : nDofRequestedFromRanks_)
  {
    int foreignRankNo = nDofRequestedFromRank.first;
    int nFromRank = nDofRequestedFromRank.second;

    if (nFromRank != 0)
    {
      requestedDofsGlobalPetscSendBufferValues[i].resize(nFromRank*nComponentsDirichletBc, 0.0);

      for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
      {
        for (int requestedDofIndex = 0; requestedDofIndex < nFromRank; requestedDofIndex++)
        {
          global_no_t requestedDofNoGlobalPetsc = requestedDofsGlobalPetsc_[i][requestedDofIndex];
          bool isLocal;
          dof_no_t requestedDofNoLocal = this->meshPartition_->getDofNoLocal(requestedDofNoGlobalPetsc, isLocal);

          requestedDofsGlobalPetscSendBufferValues[i][componentNo*nFromRank + requestedDofIndex] = this->boundaryConditionValues_[componentNo][requestedDofNoLocal];
        }
      }

      MPI_Request sendRequestValues;
      MPIUtility::handleReturnValue(MPI_Isend(requestedDofsGlobalPetscSendBufferValues[i].data(), nFromRank*nComponentsDirichletBc, MPI_DOUBLE, foreignRankNo, 1,
                                              this->meshPartition_->mpiCommunicator(), &sendRequestValues), "MPI_Isend");
      sendRequests.push_back(sendRequestValues);
    }
    i++;
  }

  // receive values
  std::vector<std::vector<double>> requestedDofsGlobalPetscReceiveBufferValues(requestDofsFromRanks_.size());

  i = 0;
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
  {
    int foreignRankNo = requestDofsFromRank.first;
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();

    requestedDofsGlobalPetscReceiveBufferValues[i].resize(nRequestedDofs*nComponentsDirichletBc);

    MPI_Request receiveRequestValues;
    MPIUtility::handleReturnValue(MPI_Irecv(requestedDofsGlobalPetscReceiveBufferValues[i].data(), nRequestedDofs*nComponentsDirichletBc, MPI_DOUBLE, foreignRankNo, 1,
                                            this->meshPartition_->mpiCommunicator(), &receiveRequestValues), "MPI_Irecv");
    receiveRequests.push_back(receiveRequestValues);
    i++;
  }

  VLOG(1) << "MPI_Waitall (3), " << sendRequests.size() << "+" << receiveRequests.size() << " requests, " << *this->meshPartition_->rankSubset();

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  // copy received dofs to new vector
  // and copy received values
  i = 0;
  for (const std::pair<int,DofsRequest> &requestDofsFromRank : requestDofsFromRanks_)
  {
    int nRequestedDofs = requestDofsFromRank.second.dofNosLocal.size();

    //const std::vector<global_no_t> &ghostDofNosGlobalPetsc = requestDofsFromRank.second.dofNosGlobalPetsc;
    const std::vector<dof_no_t> &ghostDofNosLocal = requestDofsFromRank.second.dofNosLocal;

    for (int componentNo = 0; componentNo < nComponentsDirichletBc; componentNo++)
    {
      for (int j = 0; j < nRequestedDofs; j++)
      {
        double value = requestedDofsGlobalPetscReceiveBufferValues[i][componentNo*nRequestedDofs + j];
        dof_no_t dofNoLocal = ghostDofNosLocal[j];    // this value was not sent

        boundaryConditionValues_[componentNo][dofNoLocal] = value;
      }
    }
    i++;
  }
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

  LOG(DEBUG) << "createVector \"" << this->name_ << "\", global: " << vectorCombinedWithoutDirichletDofsGlobal_ << ", local: " << vectorCombinedWithoutDirichletDofsLocal_;
}
