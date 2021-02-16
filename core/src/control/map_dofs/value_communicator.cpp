#include "control/map_dofs/value_communicator.h"

#include "easylogging++.h"
#include <mpi.h>

namespace Control
{

void ValueCommunicator::
initialize(const std::map<int,std::vector<int>> remoteDofNosAtRanks, std::vector<int> &dofNos, std::shared_ptr<Partition::RankSubset> rankSubset)
{
  rankSubset_ = rankSubset;
  MPI_Comm mpiCommunicator = rankSubset_->mpiCommunicator();

  // exchange, how many values should be sent to which rank
  VLOG(1) << "rankSubset " << rankSubset_ << ", create new window";

  int nRanks = rankSubset_->size();
  int ownRankNo = rankSubset_->ownRankNo();

  // create remote accessible memory
  int nBytes = nRanks * sizeof(int);
  int displacementUnit = sizeof(int);
  MPI_Win mpiMemoryWindow;

#ifdef USE_MPI_ALLOC
  int *remoteAccessibleMemory = nullptr;     // local memory where the other ranks can write to
  MPIUtility::handleReturnValue(MPI_Win_allocate(nBytes, displacementUnit, MPI_INFO_NULL, mpiCommunicator, (void *)&remoteAccessibleMemory, &mpiMemoryWindow), "MPI_Win_allocate");

  // clear buffer
  memset(remoteAccessibleMemory, 0, nBytes);
#else

  std::vector<int> remoteAccessibleMemory(nRanks, 0);     // local memory where the other ranks can write to
  MPIUtility::handleReturnValue(MPI_Win_create((void *)remoteAccessibleMemory.data(), nBytes, displacementUnit, MPI_INFO_NULL, mpiCommunicator, &mpiMemoryWindow), "MPI_Win_create");
#endif
  
  std::vector<int> nNodesToSendToRanks(nRanks);       // temporary buffer

  // output remoteDofNosAtRanks
#ifndef NDEBUG
  for (const std::pair<int,std::vector<int>> &remoteDofNosAtRank : remoteDofNosAtRanks)
  {
    VLOG(1) << "to rank " << remoteDofNosAtRank.first << " send: ";
    for (const int &remoteDofNo : remoteDofNosAtRank.second)
    {
      VLOG(1) << " value, store remotely at " << remoteDofNo;
    }
  }
#endif

  // put number of value to send to remote rank
  for (const std::pair<int,std::vector<int>> &remoteDofNosAtRank : remoteDofNosAtRanks)
  {
    int foreignRankNo = remoteDofNosAtRank.first;
    int nNodesToSend = remoteDofNosAtRank.second.size();
    nNodesToSendToRanks[foreignRankNo] = nNodesToSend;

    int offset = ownRankNo;
    VLOG(1) << "put value " << nNodesToSendToRanks[foreignRankNo] << " to rank " << foreignRankNo << ", offset " << offset;
    VLOG(1) << "MPI_Win_lock on rank " << foreignRankNo;

    // start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
    MPIUtility::handleReturnValue(MPI_Win_lock(MPI_LOCK_SHARED, foreignRankNo, 0, mpiMemoryWindow), "MPI_Win_lock");

    MPIUtility::handleReturnValue(MPI_Put(&nNodesToSendToRanks[foreignRankNo], 1, MPI_INT, foreignRankNo, offset, 1, MPI_INT, mpiMemoryWindow), "MPI_Put");

    MPIUtility::handleReturnValue(MPI_Win_unlock(foreignRankNo, mpiMemoryWindow), "MPI_Win_unlock");
  }

  MPIUtility::handleReturnValue(MPI_Win_fence(MPI_MODE_NOSUCCEED, mpiMemoryWindow), "MPI_Win_fence");

  // parse own memory what other ranks have written there and store in nDofsToReceiveFromRanks_
  //std::vector<std::pair<int,int>> nDofsToReceiveFromRanks_;   /// (foreignRank,nNodes), number of nodes requested by and to be send to foreignRank
  nDofsToReceiveFromRanks_.clear();
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    int nDofsToReceiveFromRank = remoteAccessibleMemory[rankNo];

    if (nDofsToReceiveFromRank > 0)
    {
      nDofsToReceiveFromRanks_.push_back(std::pair<int,int>(rankNo,nDofsToReceiveFromRank));
    }
  }

  // deallocate mpi memory
  MPIUtility::handleReturnValue(MPI_Win_free(&mpiMemoryWindow), "MPI_Win_free");

  VLOG(1) << "after fence, nDofsToReceiveFromRanks_: " << nDofsToReceiveFromRanks_;

  // send remoteDofNosAtRanks
  std::vector<MPI_Request> sendRequests;
  std::vector<std::vector<int>> sendBuffer(remoteDofNosAtRanks.size());

  // loop over rank with which to communicate
  int i = 0;
  for (typename std::map<int,std::vector<int>>::const_iterator iter = remoteDofNosAtRanks.begin(); iter != remoteDofNosAtRanks.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nValuesToSendTo = iter->second.size();

    std::copy(iter->second.begin(), iter->second.end(), std::back_inserter(sendBuffer[i]));

    assert(sendBuffer[i].size() == nValuesToSendTo);

    int tag = foreignRankNo*10000+nValuesToSendTo;
    MPI_Request sendRequest;
    MPIUtility::handleReturnValue(MPI_Isend(sendBuffer[i].data(), nValuesToSendTo, MPI_INT, foreignRankNo, tag,
                                            mpiCommunicator, &sendRequest), "MPI_Isend");

    VLOG(1) << "to rank " << foreignRankNo << " send " << nValuesToSendTo << " dof nos: " << sendBuffer[i] << ", tag=" << tag;

    sendRequests.push_back(sendRequest);
  }

  // receive the local dofs nos
  std::vector<std::vector<int>> dofNosFromRanks;   //< indexing same as in nDofsToReceiveFromRanks_, the dofNos given from that rank
  dofNosFromRanks.resize(nDofsToReceiveFromRanks_.size());
  std::vector<MPI_Request> receiveRequests;
  nValuesToReceive_ = 0;

  i = 0;
  for (typename std::vector<std::pair<int,int>>::iterator iter = nDofsToReceiveFromRanks_.begin(); iter != nDofsToReceiveFromRanks_.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nFromRank = iter->second;

    if (nFromRank != 0)
    {
      int tag = ownRankNo*10000+nFromRank;
      VLOG(1) << "i=" << i << ", from rank " << foreignRankNo << " receive " << nFromRank << " requests, tag=" << tag;

      dofNosFromRanks[i].resize(nFromRank);
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(dofNosFromRanks[i].data(), nFromRank, MPI_INT, foreignRankNo, tag,
                                              mpiCommunicator, &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);
      nValuesToReceive_ += nFromRank;
    }
  }

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  sendRequests.clear();
  receiveRequests.clear();

  // here,
  // * nDofsToReceiveFromRanks_ contains the number of dofs to receive from foreign ranks
  // * dofNosFromRanks contains for every foreign rank the dof nos where to store the dof values that will come from that rank

  // set values in dofNos
  for (std::vector<int> dofNosFromRank : dofNosFromRanks)
  {
    std::vector<int> outputVectorIndices;
    for (int dofNoFromRank: dofNosFromRank)
    {
      dofNos.push_back(dofNoFromRank);
    }
  }
}

//! communicate the values to the ranks
void ValueCommunicator::
communicate(const std::map<int,std::vector<double>> valuesToSendToRanks, std::vector<double> &receivedValues)
{
  MPI_Comm mpiCommunicator = rankSubset_->mpiCommunicator();

  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;

  // iterate over rank to send values to
  for (const std::pair<int,std::vector<double>> &valuesToSendToRank : valuesToSendToRanks)
  {
    int foreignRankNo = valuesToSendToRank.first;
    int nToRank = valuesToSendToRank.second.size();

    if (nToRank != 0)
    {
      MPI_Request sendRequest;
      MPIUtility::handleReturnValue(MPI_Isend(valuesToSendToRank.second.data(), nToRank, MPI_DOUBLE, foreignRankNo, 0,
                                              mpiCommunicator, &sendRequest), "MPI_Isend");
      sendRequests.push_back(sendRequest);
    }
  }

  // receive values
  receivedValues.resize(nValuesToReceive_);
  int receivedValuesIndex = 0;
  for (typename std::vector<std::pair<int,int>>::iterator iter = nDofsToReceiveFromRanks_.begin(); iter != nDofsToReceiveFromRanks_.end(); iter++)
  {
    int foreignRankNo = iter->first;
    int nFromRank = iter->second;

    if (nFromRank != 0)
    {
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(receivedValues.data()+receivedValuesIndex, nFromRank, MPI_DOUBLE, foreignRankNo, 0,
                                              mpiCommunicator, &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);

      receivedValuesIndex += nFromRank;
    }
  }

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  sendRequests.clear();
  receiveRequests.clear();
}

}  // namespace Control
