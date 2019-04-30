#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/vector_operators.h"
#include "control/types.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,int nComponents>
void DirichletBoundaryConditions<FunctionSpaceType,nComponents>::
initializeGhostElements()
{
  LOG(TRACE) << "initializeGhostElements";

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // determine own ghost elements that can be send to other ranks
  // loop over elements that have nodes with prescribed boundary conditions, only for those the integral term is non-zero
  for (typename std::vector<typename DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::ElementWithNodes>::const_iterator iter = this->boundaryConditionElements_.cbegin();
       iter != this->boundaryConditionElements_.cend(); iter++)
  {
    element_no_t elementNoLocal = iter->elementNoLocal;

    // get the dof nos of the current element
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

    typename std::vector<std::pair<int,std::array<double,nComponents>>>::const_iterator boundaryConditionDofIter = iter->elementalDofIndex.begin();

    VLOG(1) << "element " << elementNoLocal << ", elemental dofs: " << iter->elementalDofIndex;

    std::map<int,std::vector<global_no_t>> nonBoundaryConditionDofsOfRankGlobalPetsc;
    std::vector<global_no_t> boundaryConditionDofsGlobalPetsc;
    std::vector<ValueType> boundaryConditionValues;

    // loop over all dofs of the element
    for (int elementalDofIndex = 0; elementalDofIndex < nDofsPerElement; elementalDofIndex++)
    {
      dof_no_t boundaryConditionDofNoLocal = dofNosLocal[elementalDofIndex];
      node_no_t nodeNoLocal = boundaryConditionDofNoLocal / nDofsPerNode;
      int nodalDofIndex = boundaryConditionDofNoLocal % nDofsPerNode;

      int nextBoundaryConditionDofElementalDofIndex = -1;
      if (boundaryConditionDofIter != iter->elementalDofIndex.end())
        nextBoundaryConditionDofElementalDofIndex = boundaryConditionDofIter->first;

      VLOG(1) << " element " << elementNoLocal << ", dofIndex " << elementalDofIndex << ", dof " << boundaryConditionDofNoLocal
        << ", nextBoundaryConditionDofElementalDofIndex: " << nextBoundaryConditionDofElementalDofIndex;

      if (nextBoundaryConditionDofElementalDofIndex == elementalDofIndex)
      {
        // current dof has a prescribed boundary condition value

        global_no_t nodeNoGlobalPetsc = this->functionSpace_->meshPartition()->getNodeNoGlobalPetsc(nodeNoLocal);
        global_no_t dofNoGlobalPetsc = nodeNoGlobalPetsc*nDofsPerNode + nodalDofIndex;
        boundaryConditionDofsGlobalPetsc.push_back(dofNoGlobalPetsc);
        boundaryConditionValues.push_back(boundaryConditionDofIter->second);

        VLOG(1) << "   dof has prescribed value " << boundaryConditionDofIter->second << ", dofNoGlobalPetsc: " << dofNoGlobalPetsc;

        boundaryConditionDofIter++;
      }
      else
      {
        VLOG(1) << "    dof has no prescribed value";

        // current dof has no prescribed boundary condition value
        int neighbourRankNo = 0;

        // if the current dof is a ghost dof on the own rank
        if (!this->functionSpace_->meshPartition()->isNonGhost(nodeNoLocal, neighbourRankNo))
        {
          // current dof is non-BC, ghost
          // current dof is owned by rank neighbourRankNo
          // determine global Petsc number
          global_no_t nodeNoGlobalPetsc = this->functionSpace_->meshPartition()->getNodeNoGlobalPetsc(nodeNoLocal);
          global_no_t dofNoGlobalPetsc = nodeNoGlobalPetsc*nDofsPerNode + nodalDofIndex;

          VLOG(1) << "    dof is ghost, dofNoGlobalPetsc: " << dofNoGlobalPetsc << ", neighbourRankNo: " << neighbourRankNo;

          nonBoundaryConditionDofsOfRankGlobalPetsc[neighbourRankNo].push_back(dofNoGlobalPetsc);
        }
      }
    }

    for (const std::pair<int,std::vector<global_no_t>> &nonBoundaryConditionDofs : nonBoundaryConditionDofsOfRankGlobalPetsc)
    {
      // create new ghost element
      GhostElement ghostElement;
      ghostElement.nonBoundaryConditionDofsOfRankGlobalPetsc.assign(nonBoundaryConditionDofs.second.begin(), nonBoundaryConditionDofs.second.end());
      ghostElement.boundaryConditionDofsGlobalPetsc.assign(boundaryConditionDofsGlobalPetsc.begin(), boundaryConditionDofsGlobalPetsc.end());
      ghostElement.boundaryConditionValues.assign(boundaryConditionValues.begin(), boundaryConditionValues.end());

      int rankNo = nonBoundaryConditionDofs.first;

      // store ghost element to foreignGhostElements_
      foreignGhostElements_[rankNo].push_back(ghostElement);
    }
  }

  VLOG(1) << "determined foreignGhostElements_: ";

  for (typename std::map<int,std::vector<GhostElement>>::iterator iter = foreignGhostElements_.begin(); iter != foreignGhostElements_.end(); iter++)
  {
    VLOG(1) << "  rank " << iter->first << " has " << iter->second.size() << " ghost elements";
    for (int i = 0; i < iter->second.size(); i++)
    {
      VLOG(1) << ", non-BC: " << iter->second[i].nonBoundaryConditionDofsOfRankGlobalPetsc
        << ", BC: " << iter->second[i].boundaryConditionDofsGlobalPetsc
        << ", values: " << iter->second[i].boundaryConditionValues;
    }
  }

  // send and receive ghost elements
  int nRanks = this->functionSpace_->meshPartition()->nRanks();
  int ownRankNo = this->functionSpace_->meshPartition()->ownRankNo();
  MPI_Comm communicator = this->functionSpace_->meshPartition()->mpiCommunicator();

  // exchange number of ghost elements to send/receive, open a window for other processes to write into how many ghost elements they will send
  // open window for MPI RMA
  /*void *remoteAccessibleMemory = nullptr;
  MPI_Win mpiMemoryWindow;
  MPIUtility::handleReturnValue(MPI_Win_allocate(nRanks*sizeof(int), sizeof(int), MPI_INFO_NULL, communicator, &remoteAccessibleMemory, &mpiMemoryWindow), "MPI_Win_allocate");

  // set to 0
  memset(remoteAccessibleMemory, 0, sizeof(int)*nRanks);
*/


  LOG(DEBUG) << "rankSubset " << *this->functionSpace_->meshPartition()->rankSubset() << ", create new window";
  std::vector<int> remoteAccessibleMemory(nRanks, 0);
  int nBytes = nRanks*sizeof(int);
  int displacementUnit = sizeof(int);
  MPI_Win mpiMemoryWindow;
  MPIUtility::handleReturnValue(MPI_Win_create((void *)remoteAccessibleMemory.data(), nBytes, displacementUnit, MPI_INFO_NULL, communicator, &mpiMemoryWindow), "MPI_Win_create");

  std::vector<int> localMemory(nRanks);

  // put number of ghost elements to the corresponding processes
  for (const std::pair<int,std::vector<GhostElement>> &ghostElement : foreignGhostElements_)
  {
    int foreignRankNo = ghostElement.first;
    int nGhostElements = ghostElement.second.size();
    localMemory[foreignRankNo] = nGhostElements;

    VLOG(1) << "put value " << localMemory[foreignRankNo] << " to rank " << foreignRankNo << ", offset " << ownRankNo;

    // start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
    MPIUtility::handleReturnValue(MPI_Win_lock(MPI_LOCK_SHARED, foreignRankNo, 0, mpiMemoryWindow), "MPI_Win_lock");

    MPIUtility::handleReturnValue(MPI_Put(&localMemory[foreignRankNo], 1, MPI_INT, foreignRankNo, ownRankNo, 1, MPI_INT, mpiMemoryWindow), "MPI_Put");

    MPIUtility::handleReturnValue(MPI_Win_unlock(foreignRankNo, mpiMemoryWindow), "MPI_Win_unlock");

  }

  MPIUtility::handleReturnValue(MPI_Win_fence(MPI_MODE_NOSUCCEED, mpiMemoryWindow), "MPI_Win_fence");

  std::vector<std::pair<int,int>> nElementsFromRanks;   /// (foreignRank,nElements), number of elements to receive from foreignRank
  for (int i = 0; i < nRanks; i++)
  {
    VLOG(1) << " rank " << i << " nGhostElements: " << remoteAccessibleMemory[i];
    if (remoteAccessibleMemory[i] > 0)
    {
      nElementsFromRanks.push_back(std::pair<int,int>(i,remoteAccessibleMemory[i]));
    }
  }

  // deallocate mpi memory
  MPIUtility::handleReturnValue(MPI_Win_free(&mpiMemoryWindow), "MPI_Win_free");

  VLOG(1) << "after fence, nElementsFromRanks: " << nElementsFromRanks;

  // send lengths of arrays in ghost elements
  //int *sendBuffer[foreignGhostElements_.size()];
  std::vector<std::vector<int>> sendBuffer(foreignGhostElements_.size());
  std::vector<MPI_Request> sendRequests;
  int i = 0;
  for (typename std::map<int,std::vector<GhostElement>>::const_iterator iter = foreignGhostElements_.cbegin(); iter != foreignGhostElements_.cend(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nGhostElements = iter->second.size();

    if (nGhostElements != 0)
    {
      //sendBuffer[i] = new int [nGhostElements*2];  // for every ghostElements the sizes of nonBoundaryConditionDofsOfRankGlobalPetsc and boundaryConditionDofsGlobalPetsc
      sendBuffer[i].resize(nGhostElements*2);  // for every ghostElements the sizes of nonBoundaryConditionDofsOfRankGlobalPetsc and boundaryConditionDofsGlobalPetsc

      int j = 0;
      for (typename std::vector<GhostElement>::const_iterator ghostElementIter = iter->second.begin(); ghostElementIter != iter->second.end(); ghostElementIter++, j++)
      {
        sendBuffer[i][2*j + 0] = ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc.size();
        sendBuffer[i][2*j + 1] = ghostElementIter->boundaryConditionDofsGlobalPetsc.size();
      }

      VLOG(1) << "send to foreignRank " << foreignRankNo << " sizes ";
      for (int k = 0; k < nGhostElements*2; k++)
      {
        VLOG(1) << "  " << sendBuffer[i][k];
      }

      assert(sendBuffer[i].size() == nGhostElements*2);
      MPI_Request sendRequest;
      MPIUtility::handleReturnValue(MPI_Isend(sendBuffer[i].data(), nGhostElements*2, MPI_INT, foreignRankNo, 0, communicator, &sendRequest), "MPI_Isend");
      sendRequests.push_back(sendRequest);
    }
  }

  // receive lengths of arrays
  std::vector<std::vector<int>> receiveBuffer(nElementsFromRanks.size());   // for every ghostElements the sizes of nonBoundaryConditionDofsOfRankGlobalPetsc and boundaryConditionDofsGlobalPetsc
  std::vector<MPI_Request> receiveRequests;

  // receive data
  i = 0;
  for (std::vector<std::pair<int,int>>::iterator nElementsFromRanksIter = nElementsFromRanks.begin(); nElementsFromRanksIter != nElementsFromRanks.end(); nElementsFromRanksIter++, i++)
  {
    int foreignRankNo = nElementsFromRanksIter->first;
    int nGhostElementsFromRank = nElementsFromRanksIter->second;

    if (nGhostElementsFromRank != 0)
    {
      receiveBuffer[i].resize(nGhostElementsFromRank*2);
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(receiveBuffer[i].data(), nGhostElementsFromRank*2, MPI_INT, foreignRankNo, 0, communicator, &receiveRequest), "MPI_Irecv");
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
    i = 0;
    for (std::vector<std::pair<int,int>>::iterator nElementsFromRanksIter = nElementsFromRanks.begin(); nElementsFromRanksIter != nElementsFromRanks.end(); nElementsFromRanksIter++, i++)
    {
      int foreignRank = nElementsFromRanksIter->first;
      int nGhostElementsFromRank = nElementsFromRanksIter->second;

      if (nGhostElementsFromRank != 0)
      {
        for (int j = 0; j < nGhostElementsFromRank; j++)
        {
          VLOG(1) << "  foreignRank " << foreignRank << ", ghost element " << j << ", number nonBC dofs: " << receiveBuffer[i][2*j+0];
          VLOG(1) << "               number BC dofs: " << receiveBuffer[i][2*j+1];
        }
      }
    }
  }

  // send actual ghost elements
  std::vector<MPI_Request> sendRequestsValues;
  std::vector<std::vector<long>> sendBufferGhostElements(foreignGhostElements_.size());
  std::vector<std::vector<double>> sendBufferValues(foreignGhostElements_.size());
  i = 0;
  for (typename std::map<int,std::vector<GhostElement>>::const_iterator iter = foreignGhostElements_.cbegin(); iter != foreignGhostElements_.cend(); iter++, i++)
  {
    int foreignRankNo = iter->first;

    // determine send buffer size
    int sendBufferSize = 0;
    int sendBufferValuesSize = 0;
    for (typename std::vector<GhostElement>::const_iterator ghostElementIter = iter->second.cbegin(); ghostElementIter != iter->second.cend(); ghostElementIter++)
    {
      sendBufferSize += ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc.size();
      sendBufferSize += ghostElementIter->boundaryConditionDofsGlobalPetsc.size();
      sendBufferValuesSize += ghostElementIter->boundaryConditionDofsGlobalPetsc.size()*nComponents;
    }

    if (sendBufferSize != 0)
    {
      sendBufferGhostElements[i].resize(sendBufferSize);
      sendBufferValues[i].resize(sendBufferValuesSize);

      VLOG(1) << " send to " << foreignRankNo << ", len buffer for ghost elements: " << sendBufferSize << ", sendBufferValuesSize: " << sendBufferValuesSize;

      // fill send buffer
      int j1 = 0, j2 = 0;
      for (typename std::vector<GhostElement>::const_iterator ghostElementIter = iter->second.cbegin(); ghostElementIter != iter->second.cend(); ghostElementIter++)
      {
        for (int k = 0; k < ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc.size(); k++)
        {
          sendBufferGhostElements[i][j1++] = (long)ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc[k];
          VLOG(1) << "store ghost value (non-BC) to send buffer at j=" << j1-1 << ", value: " << ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc[k];
        }

        for (int k = 0; k < ghostElementIter->boundaryConditionDofsGlobalPetsc.size(); k++)
        {
          sendBufferGhostElements[i][j1++] = (long)ghostElementIter->boundaryConditionDofsGlobalPetsc[k];
          VLOG(1) << "store ghost value (BC) to send buffer at j=" << j1-1 << ", value: " << ghostElementIter->boundaryConditionDofsGlobalPetsc[k];
          VLOG(1) << "store value to send buffer starting at j=" << j2 << " (" << nComponents << " components), "
            << "values: " << ghostElementIter->boundaryConditionValues[k];
          for (int componentNo = 0; componentNo < nComponents; componentNo++)
          {
            sendBufferValues[i][j2++] = (double)ghostElementIter->boundaryConditionValues[k][componentNo];
          }
        }
      }

      // start send
      MPI_Request sendRequest, sendRequestValues;
      assert(sendBufferGhostElements[i].size() == sendBufferSize);
      MPIUtility::handleReturnValue(MPI_Isend(sendBufferGhostElements[i].data(), sendBufferSize, MPI_LONG, foreignRankNo, 0, communicator, &sendRequest), "MPI_Isend");
      sendRequests.push_back(sendRequest);

      assert(sendBufferValues[i].size() == sendBufferValuesSize);
      MPIUtility::handleReturnValue(MPI_Isend(sendBufferValues[i].data(), sendBufferValuesSize, MPI_DOUBLE, foreignRankNo, 0, communicator, &sendRequestValues), "MPI_Isend");
      sendRequestsValues.push_back(sendRequestValues);
    }
  }

  // receive ghost elements
  std::vector<std::vector<long>> receiveBufferGhostElements(nElementsFromRanks.size());
  std::vector<std::vector<double>> receiveBufferValues(nElementsFromRanks.size());
  std::vector<MPI_Request> receiveRequestsValues;

  // receive data
  i = 0;
  for (std::vector<std::pair<int,int>>::iterator nElementsFromRanksIter = nElementsFromRanks.begin(); nElementsFromRanksIter != nElementsFromRanks.end(); nElementsFromRanksIter++, i++)
  {
    int foreignRankNo = nElementsFromRanksIter->first;
    int nGhostElementsFromRank = nElementsFromRanksIter->second;

    // determine sizes of receive buffers
    int receiveBufferSize = 0;
    int receiveBufferValuesSize = 0;
    for (int ghostElementIndex = 0; ghostElementIndex < nGhostElementsFromRank; ghostElementIndex++)
    {
      receiveBufferSize += receiveBuffer[i][2*ghostElementIndex+0];     // number of nonBC dofs in the ghostElement on that rank
      receiveBufferSize += receiveBuffer[i][2*ghostElementIndex+1];     // number of BC dofs in the ghostElement on that rank
      receiveBufferValuesSize += receiveBuffer[i][2*ghostElementIndex+1]*nComponents;
    }

    LOG(DEBUG) << "receive from " << foreignRankNo << ", size: " << receiveBufferValuesSize;

    if (receiveBufferSize != 0)
    {
      receiveBufferGhostElements[i].resize(receiveBufferSize);
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(receiveBufferGhostElements[i].data(), receiveBufferSize, MPI_LONG, foreignRankNo, 0, communicator, &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);

      receiveBufferValues[i].resize(receiveBufferValuesSize);
      MPI_Request receiveRequestValues;
      MPIUtility::handleReturnValue(MPI_Irecv(receiveBufferValues[i].data(), receiveBufferValuesSize, MPI_DOUBLE, foreignRankNo, 0, communicator, &receiveRequestValues), "MPI_Irecv");
      receiveRequestsValues.push_back(receiveRequestValues);
    }
  }

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!sendRequestsValues.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequestsValues.size(), sendRequestsValues.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequestsValues.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequestsValues.size(), receiveRequestsValues.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");


  // store received values in ownGhostElements_
  i = 0;
  for (std::vector<std::pair<int,int>>::iterator nElementsFromRanksIter = nElementsFromRanks.begin(); nElementsFromRanksIter != nElementsFromRanks.end(); nElementsFromRanksIter++, i++)
  {
    int nGhostElementsFromRank = nElementsFromRanksIter->second;

    int j = 0;
    int j2 = 0;
    for (int ghostElementIndex = 0; ghostElementIndex < nGhostElementsFromRank; ghostElementIndex++)
    {
      int sizeNonBoundaryConditionDofsOfRankGlobalPetsc = receiveBuffer[i][2*ghostElementIndex+0];
      int sizeBoundaryConditionDofsGlobalPetsc = receiveBuffer[i][2*ghostElementIndex+1];    // number of values of type ValueType (vectors of nComponents)

      assert(sizeof(global_no_t) == sizeof(long)); // global_no_t = std::size_t == long should be 8 bytes on 64-bit linux
      GhostElement ghostElement;
      //ghostElement.nonBoundaryConditionDofsOfRankGlobalPetsc.assign(receiveBufferGhostElements[i].begin()+j, receiveBufferGhostElements[i].begin()+j+sizeNonBoundaryConditionDofsOfRankGlobalPetsc);
      assert(receiveBufferGhostElements[i].size() >= j+sizeNonBoundaryConditionDofsOfRankGlobalPetsc);
      ghostElement.nonBoundaryConditionDofsOfRankGlobalPetsc.assign(receiveBufferGhostElements[i].begin()+j, receiveBufferGhostElements[i].begin()+j+sizeNonBoundaryConditionDofsOfRankGlobalPetsc);
      j += sizeNonBoundaryConditionDofsOfRankGlobalPetsc;

      assert(receiveBufferGhostElements[i].size() >= j+sizeBoundaryConditionDofsGlobalPetsc);
      ghostElement.boundaryConditionDofsGlobalPetsc.assign(receiveBufferGhostElements[i].begin()+j, receiveBufferGhostElements[i].begin()+j+sizeBoundaryConditionDofsGlobalPetsc);
      j += sizeBoundaryConditionDofsGlobalPetsc;

      //LOG(DEBUG) << "receiveBufferValues[i=" << i << "] size: " << receiveBufferValues[i].size() << "(" << receiveBufferValues << ")" << ", j2=" << j2
      //  << ", ghostElementIndex: " << ghostElementIndex << ", sizeBoundaryConditionDofsGlobalPetsc: " << sizeBoundaryConditionDofsGlobalPetsc << " (" << j2+sizeBoundaryConditionDofsGlobalPetsc << ")";
      assert(receiveBufferValues[i].size() >= j2+sizeBoundaryConditionDofsGlobalPetsc);

      ghostElement.boundaryConditionValues.resize(sizeBoundaryConditionDofsGlobalPetsc);
      for (int valueIndex = 0; valueIndex < sizeBoundaryConditionDofsGlobalPetsc; valueIndex++)
      {
        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          ghostElement.boundaryConditionValues[valueIndex][componentNo] = receiveBufferValues[i][j2 + valueIndex*nComponents + componentNo];
        }
      }
      j2 += sizeBoundaryConditionDofsGlobalPetsc*nComponents;

      ownGhostElements_.push_back(ghostElement);
    }
  }

  VLOG(1) << "received ownGhostElements_: ";
  for (typename std::vector<GhostElement>::iterator iter = ownGhostElements_.begin(); iter != ownGhostElements_.end(); iter++)
  {
    VLOG(1) << "  non-BC: " << iter->nonBoundaryConditionDofsOfRankGlobalPetsc << ", BC: " << iter->boundaryConditionDofsGlobalPetsc
      << ", values: " << iter->boundaryConditionValues;
  }

  /*
   *  struct GhostElement
  {
    std::vector<global_no_t> nonBoundaryConditionDofsOfRankGlobalPetsc;    ///< the non-BC dofs of this element, as global petsc no. that are owned by the rank with no neighbouringRankNo
    std::vector<global_no_t> boundaryConditionDofsGlobalPetsc;      ///< the Dirichlet BC dofs of this element
    std::vector<double> boundaryConditionValues;   ///< the prescribed value, corresponding to boundaryConditionDofsGlobalPetsc
  };
  std::map<int,std::vector<GhostElement>> foreignGhostElements_;   ///< ghost elements that are normal elements on this rank, key is the rankNo of the rank to send them to
  std::vector<GhostElement> ownGhostElements_;   ///< the ghost elements for this rank
  */
}

// set the boundary conditions to system matrix, i.e. zero rows and columns of Dirichlet BC dofs and set diagonal to 1
template<typename FunctionSpaceType,int nComponents>
void DirichletBoundaryConditions<FunctionSpaceType,nComponents>::
applyInSystemMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix,
                    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> boundaryConditionsRightHandSideSummand)
{
  LOG(TRACE) << "DirichletBoundaryConditionsBase::applyInSystemMatrix";
  VLOG(1) << "boundaryConditionsRightHandSideSummand: " << *boundaryConditionsRightHandSideSummand;

  // boundary conditions for local non-ghost dofs are stored in the following member variables:
  // std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos_;        ///< vector of all local (non-ghost) boundary condition dofs
  // std::vector<ValueType> boundaryConditionValues_;               ///< vector of the local prescribed values, related to boundaryConditionDofLocalNos_

  // with typedef std::array<double,nComponents> ValueType;   ///< the type of value of one boundary condition

  // determine actions to be executed later, this is such that there are no duplicates
  // One action is the following:
  //   1. get matrix entry M_{row,col}
  //   2. update rhs rhs_{row} -= M_{row,col}*BC_col
  // The row and column indices are stored in global PETSc ordering.
  std::map<global_no_t, std::pair<ValueType, std::set<global_no_t>>> action; // map[columnNoGlobalPetsc] = <bc value, <rowNosGlobalPetsc>>

  // save matrix entries to use them later to adjust the rhs entries
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // loop over elements that have nodes with prescribed boundary conditions, only for those the integral term is non-zero
  for (typename std::vector<typename DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::ElementWithNodes>::const_iterator iter = this->boundaryConditionElements_.cbegin();
       iter != this->boundaryConditionElements_.cend(); iter++)
  {
    // get dofs of element (with and without ghosts), at least some of the dofs of this element are prescribed Dirchlet boundary condition values
    element_no_t elementNoLocal = iter->elementNoLocal;
    std::vector<dof_no_t> rowDofNosLocalWithoutGhosts;
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = this->functionSpace_->getElementDofNosLocal(elementNoLocal);
    this->functionSpace_->getElementDofNosLocalWithoutGhosts(elementNoLocal, rowDofNosLocalWithoutGhosts);

    // get row dofs global petsc nos
    std::vector<PetscInt> rowDofNosGlobalPetsc;
    this->functionSpace_->meshPartition()->getDofNoGlobalPetsc(rowDofNosLocalWithoutGhosts, rowDofNosGlobalPetsc);


    VLOG(1) << "element " << elementNoLocal << ", dofNosLocal: " << dofNosLocal << ", rowDofNosLocalWithoutGhosts: " << rowDofNosLocalWithoutGhosts;

    // loop over dofs of element that have prescribed Dirichlet boundary condition values
    for (typename std::vector<std::pair<int,ValueType>>::const_iterator columnDofsIter = iter->elementalDofIndex.begin(); columnDofsIter != iter->elementalDofIndex.end(); columnDofsIter++)
    {
      int elementalDofIndexColumn = columnDofsIter->first;
      ValueType boundaryConditionValue = columnDofsIter->second;

      dof_no_t boundaryConditionColumnDofNoLocal = dofNosLocal[elementalDofIndexColumn];

      // store the boundary condition value to action
      global_no_t boundaryConditionColumnDofNoGlobal = this->functionSpace_->meshPartition()->getDofNoGlobalPetsc(boundaryConditionColumnDofNoLocal);
      action[boundaryConditionColumnDofNoGlobal].first = boundaryConditionValue;
      action[boundaryConditionColumnDofNoGlobal].second.insert(rowDofNosGlobalPetsc.begin(), rowDofNosGlobalPetsc.end());

      // do only store action, do not perform yet, all collected actions are duplicate-cleared (automatically, because of map) and executed at the end of this method

      // commented out code would perform action
      /*
      VLOG(1) << "  dof " << boundaryConditionColumnDofNoLocal << ", BC value: " << boundaryConditionValue;
      PetscInt columnNo = boundaryConditionColumnDofNoLocal;

      // get matrix entries that correspond to column boundaryConditionColumnDofNoLocal
      std::vector<double> values(rowDofNosLocalWithoutGhosts.size());
      systemMatrix->getValues(rowDofNosLocalWithoutGhosts.size(), rowDofNosLocalWithoutGhosts.data(), 1, &columnNo, values.data());

      // scale values with -boundaryConditionValue
      for (double &v : values)
      {
        v *= -boundaryConditionValue;
      }

      VLOG(1) << "rhs set values at " << rowDofNosLocalWithoutGhosts << ", values: " << values;

      // subtract values*boundaryConditionValue from boundaryConditionsRightHandSideSummand
      boundaryConditionsRightHandSideSummand->setValues(rowDofNosLocalWithoutGhosts, values, ADD_VALUES);*/
    }
  }

  VLOG(1) << ownGhostElements_.size() << " ghost elements";

  // loop over ghost elements
  for (typename std::vector<GhostElement>::iterator ghostElementIter = ownGhostElements_.begin(); ghostElementIter != ownGhostElements_.end(); ghostElementIter++)
  {
    std::vector<PetscInt> rowDofsGlobal;
    rowDofsGlobal.assign(ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc.begin(), ghostElementIter->nonBoundaryConditionDofsOfRankGlobalPetsc.end());

    // transform the global nos to local nos
    std::vector<PetscInt> rowDofsLocal(rowDofsGlobal.size());
    std::transform(rowDofsGlobal.begin(), rowDofsGlobal.end(), rowDofsLocal.begin(), [this](global_no_t nodeNoGlobalPetsc)
    {
      return this->functionSpace_->meshPartition()->getNodeNoLocal(nodeNoGlobalPetsc);
    });

    VLOG(1) << "rowDofs (global Petsc): " << rowDofsGlobal << ", local: " << rowDofsLocal;

    // loop over dofs that are owned by this rank, these are the row dofs to consider
    int i = 0;
    for (std::vector<global_no_t>::iterator columnDofIter = ghostElementIter->boundaryConditionDofsGlobalPetsc.begin(); columnDofIter != ghostElementIter->boundaryConditionDofsGlobalPetsc.end(); columnDofIter++, i++)
    {
      PetscInt columnDofNoGlobalPetsc = *columnDofIter;
      ValueType boundaryConditionValue = ghostElementIter->boundaryConditionValues[i];

      VLOG(1) << "  dof " << columnDofNoGlobalPetsc << " (global PETSc), BC value: " << boundaryConditionValue;

      VLOG(1) << "systemMatrix->getValuesGlobalPetscIndexing(" << rowDofsGlobal << "," << columnDofNoGlobalPetsc << ")";

      // store the boundary condition value to action
      global_no_t boundaryConditionColumnDofNoGlobal = columnDofNoGlobalPetsc;
      action[boundaryConditionColumnDofNoGlobal].first = boundaryConditionValue;
      action[boundaryConditionColumnDofNoGlobal].second.insert(rowDofsGlobal.begin(), rowDofsGlobal.end());

      // do only store action, do not perform yet, all collected actions are duplicate-cleared and executed at the end of this method

      // commented out code would perform action
/*
      std::vector<double> values(rowDofsGlobal.size());
      systemMatrix->getValuesGlobalPetscIndexing(rowDofsGlobal.size(), rowDofsGlobal.data(), 1, &columnDofNoGlobalPetsc, values.data());

      // scale values with -boundaryConditionValue
      for (double &v : values)
      {
        v *= -boundaryConditionValue;
      }

      VLOG(1) << "rhs set values at " << rowDofsLocal << ", values: " << values;


      // subtract values*boundaryConditionValue from boundaryConditionsRightHandSideSummand
      boundaryConditionsRightHandSideSummand->setValues(rowDofsLocal, values, ADD_VALUES);
      */
    }
  }

  VLOG(1) << "actions: " << action;
  VLOG(1) << "rhs summand before: " << *boundaryConditionsRightHandSideSummand;

  // debugging output
  // std::map<global_no_t, std::pair<double, std::set<global_no_t>>> action; // map[columnNoGlobalPetsc] = <bc value, <rowNosGlobalPetsc>>
#if 0
  LOG(DEBUG) << "actions: ";
  for (std::map<global_no_t, std::pair<double, std::set<global_no_t>>>::iterator actionIter = action.begin(); actionIter != action.end(); actionIter++)
  {
    PetscInt columnDofNoGlobalPetsc = actionIter->first;
    double boundaryConditionValue = actionIter->second.first;

    std::vector<PetscInt> rowDofNoGlobalPetsc(actionIter->second.second.begin(), actionIter->second.second.end());

    if (actionIter->second.second.find(16) != actionIter->second.second.end())
    {
      LOG(DEBUG) << "rhs_{16} -= M_{16," << columnDofNoGlobalPetsc << "}*" << boundaryConditionValue;
    }
  }
#endif

  std::vector<double> valuesBuffer;
  std::vector<dof_no_t> dofNosBuffer;

  LOG(DEBUG) << "dirichlet boundary conditions in system matrix: execute actions, adjust rhs summand";

  // execute actions
  for (typename std::map<global_no_t, std::pair<ValueType, std::set<global_no_t>>>::iterator actionIter = action.begin();
       actionIter != action.end(); actionIter++)
  {
    PetscInt columnDofNoGlobalPetsc = actionIter->first;
    ValueType boundaryConditionValue = actionIter->second.first;

    std::vector<PetscInt> rowDofNoGlobalPetsc(actionIter->second.second.begin(), actionIter->second.second.end());

    // transform row dofs from global petsc no to local no
    std::vector<PetscInt> rowDofNosLocal(rowDofNoGlobalPetsc.size());
    std::transform(rowDofNoGlobalPetsc.begin(), rowDofNoGlobalPetsc.end(), rowDofNosLocal.begin(), [this](global_no_t nodeNoGlobalPetsc)
    {
      return this->functionSpace_->meshPartition()->getDofNoLocal(nodeNoGlobalPetsc);
    });

    // get the values of the column from the matrix
    std::vector<ValueType> values(rowDofNoGlobalPetsc.size());
    systemMatrix->template getValuesGlobalPetscIndexing<nComponents>(rowDofNoGlobalPetsc.size(), rowDofNoGlobalPetsc.data(), 1, &columnDofNoGlobalPetsc, values);

    VLOG(1) << "system matrix, col " << columnDofNoGlobalPetsc << ", rows " << rowDofNoGlobalPetsc << ", values: " << values;

    // debugging output
#if 0
    int jj = 0;
    for (std::vector<PetscInt>::iterator rowDofNosLocalIter = rowDofNosLocal.begin(); rowDofNosLocalIter != rowDofNosLocal.end(); rowDofNosLocalIter++, jj++)
    {
      if (rowDofNoGlobalPetsc[jj] == 16)
      {
        LOG(DEBUG) << "rhs_{16} -= M_{16," << columnDofNoGlobalPetsc << "}*" << boundaryConditionValue << " = " << values[jj] << "*" << boundaryConditionValue
           << " = " <<  values[jj]*boundaryConditionValue;
        debugValue -= values[jj]*boundaryConditionValue;
      }
    }
#endif

    // scale values with -boundaryConditionValue
    for (ValueType &v : values)
    {
      for (int i = 0; i < nComponents; i++)
      {
        v[i] *= -boundaryConditionValue[i];
      }
    }

    VLOG(1) << " multiplied with BC value " << boundaryConditionValue << ": " << values;

    // subtract values*boundaryConditionValue from boundaryConditionsRightHandSideSummand
    if (nComponents == 1)
    {
      // for equations with one component, simply add values
      boundaryConditionsRightHandSideSummand->setValues(rowDofNosLocal, values, ADD_VALUES);
    }
    else
    {
      // for equations with multiple components, some components in some dofs may be set to None (NaN) which indicates that they should not be touched
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        // check if column dof has a valid Dirichlet BC, else do nothing for this component
        if (!std::isfinite(boundaryConditionValue[componentNo]))
          continue;

        valuesBuffer.clear();
        dofNosBuffer.clear();
        valuesBuffer.reserve(values.size());
        dofNosBuffer.reserve(values.size());

        // copy dofs and values for valid values to buffers
        for (int i = 0; i < values.size(); i++)
        {
          if (std::isfinite(values[i][componentNo]))
          {
            dofNosBuffer.push_back(rowDofNosLocal[i]);
            valuesBuffer.push_back(values[i][componentNo]);
          }
        }

        // perform action for component
        boundaryConditionsRightHandSideSummand->setValues(componentNo, dofNosBuffer, valuesBuffer, ADD_VALUES);
      }
    }

    VLOG(1) << " after set values at " << rowDofNosLocal << ": " << *boundaryConditionsRightHandSideSummand;
  }

  VLOG(1) << "rhs summand afterwards: " << *boundaryConditionsRightHandSideSummand;


  /*
  struct GhostElement
  {
    std::vector<global_no_> nonBoundaryConditionDofsOfRankGlobalPetsc;    // the non-BC dofs of this element, as global petsc no. that are owned by the rank with no neighbouringRankNo
    std::vector<global_no_t> boundaryConditionDofsGlobalPetsc;      // the Dirichlet BC dofs of this element
  };
  std::map<int,std::vector<GhostElement>> foreignGhostElements_;   ///< ghost elements that are normal elements on this rank, key is the rankNo of the rank to send them to
  std::vector<GhostElement> ownGhostElements_;   ///< the ghost elements for this rank
  */


  // zero entries in stiffness matrix that correspond to dirichlet dofs

  //systemMatrix->assembly(MAT_FLUSH_ASSEMBLY);

  LOG(DEBUG) << "dirichlet boundary conditions in system matrix: zero entries in system matrix";

  if (nComponents > 1)
  {
    LOG(DEBUG) << "zero columns";

    // loop over columns or boundary condition dofs
    for (typename std::map<global_no_t, std::pair<ValueType, std::set<global_no_t>>>::iterator actionIter = action.begin();
        actionIter != action.end(); actionIter++)
    {
      // get the column, this is potentially a non-local dof
      PetscInt columnDofNoGlobalPetsc = actionIter->first;
      ValueType boundaryConditionValue = actionIter->second.first;

      // get the row dofs, these are all in the local range, but the numbers are global-petsc
      std::vector<PetscInt> rowDofNoGlobalPetsc(actionIter->second.second.begin(), actionIter->second.second.end());

      // for equations with multiple components, some components in some dofs may be set to None (NaN) which indicates that they should not be touched
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        // check if column dof has a valid Dirichlet BC, else do nothing for this component
        if (!std::isfinite(boundaryConditionValue[componentNo]))
          continue;

        // set entries in nested system matrices to zero
        // Actually, the rows and columns corresponding to each BC dof have to be zeroed out. Here, we set the columns to zero,
        // later, systemMatrix->zeroRowsColumns zeros out the rows only.
        // We have to use the global indexing, because the column maybe a non-local dof, but all global columns of the local rows are stored on the own rank always.

        // zero out column in all component matrices with row component "componentIndex" and column component "componentNo"
        std::vector<double> zeros(rowDofNoGlobalPetsc.size(), 0.0);

        // loop over rows of sub matrices
        for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
        {
          // do not do anything in the diagonal submatrix
          if (componentIndex == componentNo)
          {
            continue;
          }

          int zeroComponentNo = componentIndex*nComponents + componentNo;
          VLOG(1) << "zeroComponentNo: " << zeroComponentNo << ", columnDofNoGlobalPetsc: " << columnDofNoGlobalPetsc;
          systemMatrix->setValuesGlobalPetscIndexing(zeroComponentNo, rowDofNoGlobalPetsc.size(), rowDofNoGlobalPetsc.data(), 1, &columnDofNoGlobalPetsc, zeros.data(), INSERT_VALUES);
        }
      }
    }
  }

  systemMatrix->assembly(MAT_FINAL_ASSEMBLY);


  // set values of row and column of the dofs to zero and diagonal entry to 1
  //LOG(DEBUG) << "apply dirichlet BC: zeroRowsColumns at dofs " << this->boundaryConditionNonGhostDofLocalNos_;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    VLOG(1) << "zero rowsColumns componentNo " << componentNo;
    systemMatrix->zeroRowsColumns(componentNo, this->boundaryConditionsByComponent_[componentNo].dofNosLocal.size(), this->boundaryConditionsByComponent_[componentNo].dofNosLocal.data(), 1.0);
  }

  systemMatrix->assembly(MAT_FINAL_ASSEMBLY);
  VLOG(1) << "stiffness matrix after apply Dirichlet BC: " << *systemMatrix;
}

template<typename FunctionSpaceType,int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::
applyInVector(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable)
{
  //fieldVariable->setValues(this->boundaryConditionNonGhostDofLocalNos_, this->boundaryConditionValues_);

  // set boundary condition dofs to prescribed values, only non-ghost dofs
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // std::array<std::pair<std::vector<dof_no_t>, std::vector<double>>, nComponents> boundaryConditionsByComponent_;   ///< the boundary condition data organized by component
    fieldVariable->setValues(componentNo, this->boundaryConditionsByComponent_[componentNo].dofNosLocal,
                             this->boundaryConditionsByComponent_[componentNo].values, INSERT_VALUES);

    VLOG(1) << "applied boundary conditions, component " << componentNo
      << ", local dofs: " << this->boundaryConditionsByComponent_[componentNo].dofNosLocal
      << ", values: " << this->boundaryConditionsByComponent_[componentNo].values;
  }
}

// set the boundary conditions to the right hand side
template<typename FunctionSpaceType,int nComponents>
void DirichletBoundaryConditions<FunctionSpaceType,nComponents>::
applyInRightHandSide(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSide,
                     std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> boundaryConditionsRightHandSideSummand)
{
  //LOG(TRACE) << "DirichletBoundaryConditionsBase::applyInRightHandSide";

  VLOG(1) << "applyInRightHandSide: rightHandSide=" << *rightHandSide;
  VLOG(1) << "boundaryConditionsRightHandSideSummand=" << *boundaryConditionsRightHandSideSummand;

  if (rightHandSide != boundaryConditionsRightHandSideSummand)
  {
    // set rhs += summand, where summand is the helper variable boundaryConditionsRightHandSideSummand
    PetscErrorCode ierr;
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      ierr = VecAXPY(rightHandSide->valuesGlobal(componentNo), 1, boundaryConditionsRightHandSideSummand->valuesGlobal(componentNo)); CHKERRV(ierr);
    }
  }

  this->applyInVector(rightHandSide);

  VLOG(1) << "rightHandSide after set values: " << *rightHandSide;
}


}  // namespace
