#include "partition/mesh_partition/01_mesh_partition_composite.h"

#include <iterator>

#include "easylogging++.h"
#include "utility/string_utility.h"
#include "utility/vector_operators.h"

namespace Partition
{

const double nodePositionEqualTolerance = 1e-5;

template<int D, typename BasisFunctionType>
MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
MeshPartition(const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), subFunctionSpaces_(subFunctionSpaces)
{
  nSubMeshes_ = this->subFunctionSpaces_.size();

  // initialize number of local and global elements
  this->nElementsLocal_ = 0;
  this->nElementsGlobal_ = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    // count number of elements
    this->nElementsLocal_ += subFunctionSpace->nElementsLocal();
    this->nElementsGlobal_ += subFunctionSpace->nElementsGlobal();
    LOG(DEBUG) << "sub function space \"" << subFunctionSpace->meshName() << "\" has " << subFunctionSpace->nElementsLocal() << " local, "
      << subFunctionSpace->nElementsGlobal() << " global elements.";
  }

  elementNoGlobalBegin_ = 0;
  global_no_t nElementsLocalValue = nElementsLocal_;
  MPIUtility::handleReturnValue(MPI_Exscan(&nElementsLocalValue, &elementNoGlobalBegin_, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, this->mpiCommunicator()), "MPI_Exscan");

  LOG(DEBUG) << "total elements local: " << this->nElementsLocal_ << ", global: " << this->nElementsGlobal_
    << ", elementNoGlobalBegin_: " << elementNoGlobalBegin_;


  // determine nodes that are the same on multiple meshes
  initializeSharedNodes();

  assert(removedSharedNodes_.size() == nSubMeshes_);

  // count number of shared local nodes, without ghosts and with ghosts
  nNodesSharedLocal_ = 0;
  nGhostNodesSharedLocal_ = 0;
  nRemovedNodesNonGhost_.resize(nSubMeshes_);

  // iterate over meshes
  for (int meshIndex = 0; meshIndex < nSubMeshes_; meshIndex++)
  {
    // removedSharedNodes_[meshNo][nodeNo] = <sameAsInMeshNo,nodeNoOfThatMesh>
    std::map<node_no_t,std::pair<int,node_no_t>> &sharedNodesInMesh = removedSharedNodes_[meshIndex];

    nRemovedNodesNonGhost_[meshIndex] = 0;

    // iterate over nodes of this mesh that have an equal node in other meshes
    for (std::map<node_no_t,std::pair<int,node_no_t>>::iterator iter2 = sharedNodesInMesh.begin(); iter2 != sharedNodesInMesh.end(); iter2++)
    {
      node_no_t sharedNodeNo = iter2->first;

      // if node is non-ghost
      if (sharedNodeNo < subFunctionSpaces_[meshIndex]->nNodesLocalWithoutGhosts())
      {

#ifndef NDEBUG
        // check if node positions equal
        const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
        Vec3 position0 = subFunctionSpaces_[meshIndex]->geometryField().getValue(sharedNodeNo*nDofsPerNode);
        Vec3 position1 = subFunctionSpaces_[iter2->second.first]->geometryField().getValue(iter2->second.second*nDofsPerNode);
        double distance = MathUtility::distance<3>(position0, position1);

        LOG(DEBUG) << "[meshNo " << meshIndex << ", nodeNo " << sharedNodeNo << " (non-ghost)] " << position0
          << " equals [meshNo " << iter2->second.first << ", nodeNo " << iter2->second.second << "] " << position1
          << ", distance=" << distance;

        //assert (distance < nodePositionEqualTolerance);
#endif

        nRemovedNodesNonGhost_[meshIndex]++;
        nNodesSharedLocal_++;
      }
      else
      {
        LOG(DEBUG) << "[meshNo " << meshIndex << ", nodeNo " << sharedNodeNo << " (ghost)] equals [meshNo " << iter2->second.first << ", nodeNo " << iter2->second.second << "]";

        nGhostNodesSharedLocal_++;
      }
    }
  }

  LOG(DEBUG) << "nRemovedNodesNonGhost_: " << nRemovedNodesNonGhost_;
  LOG(DEBUG) << "nNodesSharedLocal_: " << nNodesSharedLocal_ << ", nGhostNodesSharedLocal_: " << nGhostNodesSharedLocal_;

  initializeGhostNodeNos();

  createLocalDofOrderings();

  // print warning if composite mesh does not overlap

  for (int subMeshNo = 1; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    if (removedSharedNodes_[subMeshNo].empty())
    {
      LOG(WARNING) << "In composite mesh: sub mesh " << subMeshNo << " \"" << subFunctionSpaces_[subMeshNo]->meshName()
        << "\" (" << subFunctionSpaces_[subMeshNo]->nNodesGlobal() << " global nodes) does not overlap with any other submeshes.";
    }
  }

  // checking
#ifndef NDEBUG

  for (int subMeshNo = 1; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    LOG(DEBUG) << "sub mesh " << subMeshNo << " \"" << subFunctionSpaces_[subMeshNo]->meshName()
      << "\" (" << subFunctionSpaces_[subMeshNo]->nNodesGlobal() << " global nodes) shares " << removedSharedNodes_[subMeshNo].size() << " nodes";
  }

  int nRemovedNodesNonGhostTotal = 0;
  for (int i = 0; i < nRemovedNodesNonGhost_.size(); i++)
  {
    nRemovedNodesNonGhostTotal += nRemovedNodesNonGhost_[i];
  }
  assert(nNodesSharedLocal_ == nRemovedNodesNonGhostTotal);

  int nNodesLocal = 0;
  for (int i = 0; i < nRemovedNodesNonGhost_.size(); i++)
  {
    nNodesLocal += nNonDuplicateNodesWithoutGhosts_[i];
  }
  assert(nNodesLocalWithoutGhosts_ == nNodesLocal);

  // count number of ghost and non-ghost local dofs
  int nNodesLocalWithGhosts = 0;
  LOG(DEBUG) << "meshAndNodeNoLocalToNodeNoNonDuplicateLocal: " << meshAndNodeNoLocalToNodeNoNonDuplicateLocal_;
  for (const std::vector<node_no_t> &meshAndNodeNoNonDuplicateLocal : meshAndNodeNoLocalToNodeNoNonDuplicateLocal_)
  {
    for (const node_no_t nodeNoNonDuplicateLocal : meshAndNodeNoNonDuplicateLocal)
    {
      if (nodeNoNonDuplicateLocal != -1)
      {
        nNodesLocalWithGhosts++;
      }
    }
  }
  LOG(DEBUG) << "nNodesLocalWithGhosts: " << nNodesLocalWithGhosts << "=" << nNodesLocalWithGhosts_;
  assert (nNodesLocalWithGhosts == nNodesLocalWithGhosts_);

  // count number of ghost and non-ghost local dofs
  nNodesLocalWithGhosts = 0;
  for (const std::vector<node_no_t> &meshAndNodeNoNonDuplicateLocal : meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_)
  {
    for (const node_no_t nodeNoNonDuplicateLocal : meshAndNodeNoNonDuplicateLocal)
    {
      if (nodeNoNonDuplicateLocal != -1)
      {
        nNodesLocalWithGhosts++;
      }
    }
  }
  assert (nNodesLocalWithGhosts == nNodesLocalWithGhosts_);


  LOG(DEBUG) << "getString: ";
  LOG(DEBUG) << getString();
#endif

}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
initializeSharedNodes()
{
  // determine nodes that are the same on multiple meshes

  // get node positions of all meshes
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  std::vector<std::vector<Vec3>> nodePositions(nSubMeshes_);
  std::vector<std::vector<std::pair<Vec3,node_no_t>>> nodePositionsNodes(nSubMeshes_);  // the node positions with nodes for every submesh

  // iterate over submeshes and save all node positions
  int i = 0;
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    // get geometry field
    subFunctionSpace->geometryField().getValuesWithGhosts(nodePositions[i]);

    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpace->nNodesLocalWithGhosts(); nodeNoLocal++)
    {
      nodePositionsNodes[i].push_back(std::make_pair(nodePositions[i][nodeNoLocal*nDofsPerNode + 0], nodeNoLocal));
    }

    // sort according to x coordinate of node positions
    std::sort(nodePositionsNodes[i].begin(), nodePositionsNodes[i].end(), [](const std::pair<Vec3,node_no_t> &a, const std::pair<Vec3,node_no_t> &b)
    {
      return a.first[0] < b.first[0];
    });

    LOG(DEBUG) << "nodePositionsNodes of subFunctionSpace " << i << ", sorted by x: " << nodePositionsNodes[i];
    i++;
  }

  removedSharedNodes_.resize(nSubMeshes_);
  // std::vector<std::map<node_no_t,std::pair<int,node_no_t>>>

  // iterate over submeshes
  i = 0;
  for(typename std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>>::const_iterator iter = subFunctionSpaces_.cbegin();
      iter != subFunctionSpaces_.cend(); iter++, i++)
  {
    // find shared nodes in all next meshes

    LOG(DEBUG) << "find nodes of subFunctionSpace " << i << " that are the same in further meshes " << i+1 << " to " << subFunctionSpaces_.size()-1;

    if (i+1 < nSubMeshes_)
    {
      // loop over node positions of this mesh
      for (node_no_t nodeNoLocal = 0; nodeNoLocal < nodePositions[i].size(); nodeNoLocal++)
      {
        Vec3 position = nodePositions[i][nodeNoLocal];
        LOG(DEBUG) << "nodeNo " << nodeNoLocal << ", position " << position << " at x=" << position[0];

        // iterate over further submeshes
        for(int indexOtherMesh = i+1; indexOtherMesh < nSubMeshes_; indexOtherMesh++)
        {
          // find node that has closest x coordinate
          LOG(DEBUG) << "  otherMesh " << indexOtherMesh << ", find node that has the closest x coordinate to " << position[0];
          LOG(DEBUG) << "  get last node with x coordinate that is lower than " << position[0] << " by more than tolerance " << nodePositionEqualTolerance;

          // get last node with x coordinate that is lower by more than tolerance
          int k = nodePositionsNodes[indexOtherMesh].size() / 2;
          int kPrevious = -1;
          int lower = 0;
          int upper = nodePositionsNodes[indexOtherMesh].size();

          while (k != kPrevious)
          {
            Vec3 currentNodePosition = nodePositionsNodes[indexOtherMesh][k].first;
            if (currentNodePosition[0] < position[0]-nodePositionEqualTolerance)
            {
              lower = k;
            }
            else
            {
              upper = k;
            }
            kPrevious = k;
            k = (upper + lower) / 2;

            LOG(DEBUG) << "  range [" << lower << "," << upper << "] k:" << k << ", x:" << currentNodePosition[0];
          }
          LOG(DEBUG) << "  now check all node positions of otherMesh " << indexOtherMesh << " that have x=" << position[0] << " within tolerance";

          // check all node positions after k
          for (;k < nodePositionsNodes[indexOtherMesh].size(); k++)
          {
            Vec3 nodePositionOtherMesh = nodePositionsNodes[indexOtherMesh][k].first;
            node_no_t nodeNoLocalOtherMesh = nodePositionsNodes[indexOtherMesh][k].second;

            if (nodePositionOtherMesh[0] > position[0]+nodePositionEqualTolerance)
            {
              LOG(DEBUG) << "  node k: " << k << ", nodeNo: " << nodeNoLocalOtherMesh << ", position: " << nodePositionOtherMesh
                << " is over " << position[0]+nodePositionEqualTolerance << " -> break";
              break;
            }

            double distance = MathUtility::distance<3>(position, nodePositionOtherMesh);
            LOG(DEBUG) << "  node k: " << k << ", nodeNo: " << nodeNoLocalOtherMesh << ", position: " << nodePositionOtherMesh << ", distance: " << distance;

            // if the other mesh node is at the same position as the first node
            if (distance <= nodePositionEqualTolerance)
            {
              LOG(DEBUG) << "   node is shared.";
              if (removedSharedNodes_[indexOtherMesh].find(nodeNoLocal) == removedSharedNodes_[indexOtherMesh].end())
              {
                LOG(DEBUG) << "   node is not yet included in removedSharedNodes_, set (indexOtherMesh,nodeNoLocalOtherMesh) = ("
                  << i << "," << nodeNoLocal << ") at [" << indexOtherMesh << "][" << nodeNoLocalOtherMesh << "]";
                removedSharedNodes_[indexOtherMesh][nodeNoLocalOtherMesh] = std::make_pair(i, nodeNoLocal);
              }
              break;
            }
          }
        }
      }
    }
  }

  LOG(DEBUG) << "removedSharedNodes_: " << removedSharedNodes_;
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
initializeGhostNodeNos()
{
  // input values that have to be set before this initialize() method:
  // nRemovedNodesNonGhost_
  // removedSharedNodes_

  // determine the local number of non-duplicate nodes for all meshes together
  nNodesLocalWithoutGhosts_ = 0;
  nNonDuplicateNodesWithoutGhosts_.resize(nSubMeshes_);

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    int nDuplicateNodes = nRemovedNodesNonGhost_[subMeshNo];
    int nNodesLocalWithoutGhosts = subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts();
    nNonDuplicateNodesWithoutGhosts_[subMeshNo] = nNodesLocalWithoutGhosts - nDuplicateNodes;
    nNodesLocalWithoutGhosts_ += nNonDuplicateNodesWithoutGhosts_[subMeshNo];
  }

  // determine global number of duplicate-free nodes over all submeshes
  nNodesGlobal_ = 0;
  global_no_t nNodesLocalValue = nNodesLocalWithoutGhosts_;
  MPIUtility::handleReturnValue(MPI_Allreduce(&nNodesLocalValue, &nNodesGlobal_, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, this->mpiCommunicator()), "MPI_Allreduce");

  // get number of duplicate-free nodes on previous ranks
  nonDuplicateNodeNoGlobalBegin_ = 0;
  MPIUtility::handleReturnValue(MPI_Exscan(&nNodesLocalValue, &nonDuplicateNodeNoGlobalBegin_, 1, MPI_INT, MPI_SUM, this->mpiCommunicator()), "MPI_Exscan");

  // setup duplicate-free numberings (local and global) for non-ghost nodes, also store isDuplicate_ for non-ghost nodes
  meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_.resize(nSubMeshes_);
  meshAndNodeNoLocalToNodeNoNonDuplicateLocal_.resize(nSubMeshes_);
  isDuplicate_.resize(nSubMeshes_);

  global_no_t nodeNoNonDuplicateGlobal = nonDuplicateNodeNoGlobalBegin_;

  LOG(DEBUG) << "nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts_ << ", nNodesGlobal: " << nNodesGlobal_
    << ", nonDuplicateNodeNoGlobalBegin_: " << nonDuplicateNodeNoGlobalBegin_;
  LOG(DEBUG) << "now initialize numberings";

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    LOG(DEBUG) << "subMeshNo " << subMeshNo;

    int nNodesLocalWithGhosts = subFunctionSpaces_[subMeshNo]->nNodesLocalWithGhosts();
    meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo].resize(nNodesLocalWithGhosts);
    meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo].resize(nNodesLocalWithGhosts);
    isDuplicate_[subMeshNo].resize(nNodesLocalWithGhosts, false);
    nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_.reserve(nNodesLocalWithoutGhosts_);   // will contain also the entries for ghost nodes

    std::map<node_no_t,std::pair<int,node_no_t>>::iterator removedSharedNodesLocalNosIter = removedSharedNodes_[subMeshNo].begin();

    // loop over local nodes of current submesh
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      VLOG(1) << "subMeshNo " << subMeshNo << " nodeNoLocal " << nodeNoLocal;

      // if node is to be removed from the composite numbering
      if (removedSharedNodesLocalNosIter != removedSharedNodes_[subMeshNo].end())
      {
        VLOG(1) << " next removed node: " << *removedSharedNodesLocalNosIter;

        if (removedSharedNodesLocalNosIter->first == nodeNoLocal)
        {
          VLOG(1) << " -> is to be removed";

          // -1 indicates that this node is removed in the nonDuplicate numbering
          meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal] = -1;
          meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal] = -1;
          isDuplicate_[subMeshNo][nodeNoLocal] = true;

          removedSharedNodesLocalNosIter++;

          continue;
        }
      }
      VLOG(1) << " is not removed, assign new no. " << nodeNoNonDuplicateGlobal;

      // here, nodeNoLocal is a normal, non-duplicate node
      meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal] = nodeNoNonDuplicateGlobal;

      node_no_t nodeNoNonDuplicateLocal = nodeNoNonDuplicateGlobal - nonDuplicateNodeNoGlobalBegin_;
      meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal] = nodeNoNonDuplicateLocal;

      VLOG(1) << "   g: " << nodeNoNonDuplicateGlobal
        << ", nodeNoLocal " << nodeNoLocal << ", nodeNoNonDuplicateLocal: " << nodeNoNonDuplicateLocal
        << " = " << nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_.size();

      assert(nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_.size() == nodeNoNonDuplicateLocal);

      nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_.push_back(std::make_pair(subMeshNo, nodeNoLocal));

      nodeNoNonDuplicateGlobal++;
    }
  }

  VLOG(1) << "meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_ (only non-ghost nodes are set so far): " << meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_;
  VLOG(1) << "meshAndNodeNoLocalToNodeNoNonDuplicateLocal_ (only non-ghost nodes are set so far): " << meshAndNodeNoLocalToNodeNoNonDuplicateLocal_;

  // next, we communicate with neighbors to obtain duplicate-free numbering for ghost nodes

  // get vector of ghost nodes that are needed from neighbouring ranks
  /**
    * struct NodesRequest
    * {
    *   std::vector<global_no_t> nodeNosGlobalPetsc;   //< global node no
    *   std::vector<node_no_t> nodeNosLocal;            //< local node no on own rank
    * };
    *
    * std::map<int, std::vector<NodesRequest>> requestNodesFromRanks_;    //< requestNodesFromRanks_[rankNo][subMeshNo].nodeNosGlobalPetsc, for some other ranks which nodes are requested from them, for each submesh
    */

  // iterate over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    int nNodesLocalWithoutGhosts = subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts();
    int nNodesLocalWithGhosts = subFunctionSpaces_[subMeshNo]->nNodesLocalWithGhosts();

    // loop over ghost nodes of current submesh
    for (node_no_t nodeNoLocal = nNodesLocalWithoutGhosts; nodeNoLocal < nNodesLocalWithGhosts; nodeNoLocal++)
    {
      global_no_t nodeNoGlobalPetsc = subFunctionSpaces_[subMeshNo]->meshPartition()->getNodeNoGlobalPetsc(nodeNoLocal);

      // if the current node is a ghost node on the own rank (should be)
      int neighbourRankNo = 0;
      if (!subFunctionSpaces_[subMeshNo]->meshPartition()->isNonGhost(nodeNoLocal, neighbourRankNo))
      {
        // allocate space for requests for every submesh
        if (requestNodesFromRanks_.find(neighbourRankNo) == requestNodesFromRanks_.end())
        {
          requestNodesFromRanks_[neighbourRankNo].resize(nSubMeshes_);
        }

        requestNodesFromRanks_[neighbourRankNo][subMeshNo].nodeNosLocal.push_back(nodeNoLocal);
        requestNodesFromRanks_[neighbourRankNo][subMeshNo].nodeNosGlobalPetsc.push_back(nodeNoGlobalPetsc);
      }
      else
      {
        LOG(FATAL) << "ghost node not recognized as ghost node (isNonGhost is errorneous)";
      }
    }
  }

  // exchange, how many values should be sent to which rank
  VLOG(1) << "rankSubset " << *this->rankSubset() << ", create new window";

  int nRanks = this->nRanks();
  int ownRankNo = this->ownRankNo();

  // create remote accessible memory
  std::vector<int> remoteAccessibleMemory(nRanks*nSubMeshes_, 0);
  int nBytes = nRanks * nSubMeshes_ * sizeof(int);
  int displacementUnit = sizeof(int);
  MPI_Win mpiMemoryWindow;
  MPIUtility::handleReturnValue(MPI_Win_create((void *)remoteAccessibleMemory.data(), nBytes, displacementUnit, MPI_INFO_NULL, this->mpiCommunicator(), &mpiMemoryWindow), "MPI_Win_create");

  std::vector<int> localMemory(nRanks * nSubMeshes_);

  // put number of requested ghost nodes to the corresponding processes
  for (const std::pair<int,std::vector<NodesRequest>> &requestNodesFromRank : requestNodesFromRanks_)
  {
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
    {
      int foreignRankNo = requestNodesFromRank.first;
      int nRequestedNodes = requestNodesFromRank.second[subMeshNo].nodeNosLocal.size();
      localMemory[foreignRankNo*nSubMeshes_ + subMeshNo] = nRequestedNodes;

      int offset = ownRankNo*nSubMeshes_ + subMeshNo;
      VLOG(1) << "put value " << localMemory[foreignRankNo*nSubMeshes_ + subMeshNo] << " to rank " << foreignRankNo << ", offset " << offset;

      // start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
      MPIUtility::handleReturnValue(MPI_Win_lock(MPI_LOCK_SHARED, foreignRankNo, 0, mpiMemoryWindow), "MPI_Win_lock");

      MPIUtility::handleReturnValue(MPI_Put(&localMemory[foreignRankNo*nSubMeshes_ + subMeshNo], 1, MPI_INT, foreignRankNo, offset, 1, MPI_INT, mpiMemoryWindow), "MPI_Put");

      MPIUtility::handleReturnValue(MPI_Win_unlock(foreignRankNo, mpiMemoryWindow), "MPI_Win_unlock");
    }
  }

  MPIUtility::handleReturnValue(MPI_Win_fence(MPI_MODE_NOSUCCEED, mpiMemoryWindow), "MPI_Win_fence");

  // parse own memory what other ranks have written there and store in nNodesRequestedFromRanks
  std::vector<std::pair<int,std::vector<int>>> nNodesRequestedFromRanks;   /// (foreignRank,nNodes[subMeshNo]), number of nodes requested by and to be send to foreignRank
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    std::vector<int> nNodesRequestedForSubmeshes(nSubMeshes_);
    bool nodesRequestedFromThisRank = false;

    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
    {
      nNodesRequestedForSubmeshes[subMeshNo] = remoteAccessibleMemory[rankNo*nSubMeshes_ + subMeshNo];
      if (nNodesRequestedForSubmeshes[subMeshNo] > 0)
        nodesRequestedFromThisRank = true;

      VLOG(1) << " rank " << rankNo << " subMeshNo " << subMeshNo << " nRequestedNodes: " << nNodesRequestedForSubmeshes[subMeshNo];
    }

    if (nodesRequestedFromThisRank)
    {
      nNodesRequestedFromRanks.push_back(std::pair<int,std::vector<int>>(rankNo,nNodesRequestedForSubmeshes));
    }
  }

  // deallocate mpi memory
  MPIUtility::handleReturnValue(MPI_Win_free(&mpiMemoryWindow), "MPI_Win_free");

  VLOG(1) << "after fence, nNodesRequestedFromRanks: " << nNodesRequestedFromRanks;

  // exchange which nodes are requested
  std::vector<MPI_Request> sendRequests;
  std::vector<std::vector<int>> sendBuffer(requestNodesFromRanks_.size());

  int i = 0;
  for (typename std::map<int,std::vector<NodesRequest>>::iterator iter = requestNodesFromRanks_.begin(); iter != requestNodesFromRanks_.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nRequestedNodes = 0;
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      nRequestedNodes += iter->second[subMeshNo].nodeNosLocal.size();

    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      std::copy(iter->second[subMeshNo].nodeNosGlobalPetsc.begin(), iter->second[subMeshNo].nodeNosGlobalPetsc.end(), std::back_inserter(sendBuffer[i]));

    assert(sendBuffer[i].size() == nRequestedNodes);

    MPI_Request sendRequest;
    MPIUtility::handleReturnValue(MPI_Isend(sendBuffer[i].data(), nRequestedNodes, MPI_INT, foreignRankNo, 0,
                                            this->mpiCommunicator(), &sendRequest), "MPI_Isend");

    VLOG(1) << "to rank " << foreignRankNo << " send " << nRequestedNodes << " requests: " << sendBuffer[i];

    sendRequests.push_back(sendRequest);
  }

  // receive which nodes are requested
  std::vector<std::vector<int>> requestedNodesGlobalPetsc;   //< indexing same as in nNodesRequestedFromRanks_, the requested nodes from that rank
  requestedNodesGlobalPetsc.resize(nNodesRequestedFromRanks.size());
  std::vector<MPI_Request> receiveRequests;

  i = 0;
  for (typename std::vector<std::pair<int,std::vector<int>>>::iterator iter = nNodesRequestedFromRanks.begin(); iter != nNodesRequestedFromRanks.end(); iter++, i++)
  {
    int foreignRankNo = iter->first;
    int nFromRank = 0;
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      nFromRank += iter->second[subMeshNo];

    if (nFromRank != 0)
    {
      VLOG(1) << "i=" << i << ", from rank " << foreignRankNo << " receive " << nFromRank << " requests";

      requestedNodesGlobalPetsc[i].resize(nFromRank);
      MPI_Request receiveRequest;
      MPIUtility::handleReturnValue(MPI_Irecv(requestedNodesGlobalPetsc[i].data(), nFromRank, MPI_INT, foreignRankNo, 0,
                                              this->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
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
/*
  if (VLOG_IS_ON(1))
  {
    std::stringstream s;
    for (const std::pair<int,std::vector<NodesRequest>> &requestNodesFromRank : requestNodesFromRanks_)
    {
      s << "[rank " << requestNodesFromRank.first << ", nodeNosGlobalPetsc: ";
      for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      {
        s << requestNodesFromRank.second[subMeshNo].nodeNosGlobalPetsc << ", ";
      }
      s << "nodeNosLocal: ";

      for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      {
        s << requestNodesFromRank.second[subMeshNo].nodeNosLocal << ", ";
      }
      s << "], ";
    }
    VLOG(1) << "requestNodesFromRanks_: " << s.str();
    VLOG(1) << "requestedNodesGlobalPetsc: " << requestedNodesGlobalPetsc;
  }
*/
  // send nodes in nonDuplicateGlobal ordering
  i = 0;
  std::vector<std::vector<int>>    requestedNodesGlobalPetscSendBuffer(nNodesRequestedFromRanks.size());

  for (const std::pair<int,std::vector<int>> &nNodesRequestedFromRank : nNodesRequestedFromRanks)
  {
    int foreignRankNo = nNodesRequestedFromRank.first;
    int nFromRank = 0;
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      nFromRank += nNodesRequestedFromRank.second[subMeshNo];

    if (nFromRank != 0)
    {
      // set the requested global petsc nos in requestedNodesGlobalPetscSendBuffer[i]
      int requestedNodeIndex = 0;
      for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      {
        for (int j = 0; j < nNodesRequestedFromRank.second[subMeshNo]; j++, requestedNodeIndex++)
        {
          global_no_t requestedNodeNoGlobalPetsc = requestedNodesGlobalPetsc[i][requestedNodeIndex];
          bool isLocal;
          node_no_t requestedNodeNoLocal = subFunctionSpaces_[subMeshNo]->meshPartition()->getNodeNoLocal(requestedNodeNoGlobalPetsc, isLocal);
          int nodeNoNonDuplicatGlobalPetsc = meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][requestedNodeNoLocal];

          requestedNodesGlobalPetscSendBuffer[i].push_back(nodeNoNonDuplicatGlobalPetsc);
          // meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][requestedNodeNoLocal] is -1 if it is a removed node

          VLOG(1) << " send to rank " << foreignRankNo << ", subMeshNo " << subMeshNo << " requested nodeNoGlobalPetsc: " << requestedNodeNoGlobalPetsc
            << ", requestedNodeNoLocal: " << requestedNodeNoLocal << " send node " << nodeNoNonDuplicatGlobalPetsc
            << " at sendBuffer index " << requestedNodeIndex;
        }
      }
      assert(requestedNodeIndex == nFromRank);

      MPI_Request sendRequestNodes;
      MPIUtility::handleReturnValue(MPI_Isend(requestedNodesGlobalPetscSendBuffer[i].data(), nFromRank, MPI_INT, foreignRankNo, 0,
                                              this->mpiCommunicator(), &sendRequestNodes), "MPI_Isend");
      sendRequests.push_back(sendRequestNodes);
    }
    i++;
  }

  // receive nodes in nonBCGlobal ordering
  std::vector<std::vector<int>> requestedNodesGlobalPetscReceiveBuffer(requestNodesFromRanks_.size());

  i = 0;
  for (const std::pair<int,std::vector<NodesRequest>> &requestNodesFromRank : requestNodesFromRanks_)
  {
    int foreignRankNo = requestNodesFromRank.first;
    int nRequestedNodes = 0;
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      nRequestedNodes += requestNodesFromRank.second[subMeshNo].nodeNosLocal.size();

    requestedNodesGlobalPetscReceiveBuffer[i].resize(nRequestedNodes);

    MPI_Request receiveRequestNodes;
    MPIUtility::handleReturnValue(MPI_Irecv(requestedNodesGlobalPetscReceiveBuffer[i].data(), nRequestedNodes, MPI_INT, foreignRankNo, 0,
                                            this->mpiCommunicator(), &receiveRequestNodes), "MPI_Irecv");
    receiveRequests.push_back(receiveRequestNodes);
    i++;
  }

  // wait for communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  node_no_t nodeNoNonDuplicateLocal = nNodesLocalWithoutGhosts_;

  // copy received nodes to new vector
  i = 0;
  nonDuplicateGhostNodeNosGlobal_.clear();
  for (const std::pair<int,std::vector<NodesRequest>> &requestNodesFromRank : requestNodesFromRanks_)
  {
    int nRequestedNodes = 0;
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
      nRequestedNodes += requestNodesFromRank.second[subMeshNo].nodeNosLocal.size();

    int receiveBufferIndex = 0;
    for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
    {
      const std::vector<node_no_t> &ghostNodeNosLocal = requestNodesFromRank.second[subMeshNo].nodeNosLocal;
      const int nRequestedNodesOnMesh = requestNodesFromRank.second[subMeshNo].nodeNosLocal.size();

      for (int j = 0; j < nRequestedNodesOnMesh; j++, receiveBufferIndex++)
      {
        global_no_t nodeNoNonDuplicateGlobal = requestedNodesGlobalPetscReceiveBuffer[i][receiveBufferIndex];
        node_no_t nodeNoLocal = ghostNodeNosLocal[j];    // this value was not sent, it is stored in the request

        // if the received nodeNoNonDuplicateGlobal does not contain the index, but the value -1, this means the node is removed (shared with another node)
        if (nodeNoNonDuplicateGlobal == -1)
        {
          isDuplicate_[subMeshNo][nodeNoLocal] = true;
          meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal] = -1;
        }
        else
        {
          // global and local mapping
          nonDuplicateGhostNodeNosGlobal_.push_back(nodeNoNonDuplicateGlobal);
          meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal] = nodeNoNonDuplicateLocal;

          // inverse mapping
          nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_.push_back(std::make_pair(subMeshNo, nodeNoLocal));
          assert(nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_.size() == nodeNoNonDuplicateLocal);

          nodeNoNonDuplicateLocal++;
        }
        meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal] = nodeNoNonDuplicateGlobal;

        VLOG(1) << "received from rank " << requestNodesFromRank.first << ", subMeshNo " << subMeshNo << ", j: " << j << ", i: " << i
         << ", " << requestedNodesGlobalPetscReceiveBuffer[i].size() << " entries, index=" << receiveBufferIndex
         << " nodeNoLocal: " << nodeNoLocal << ", nodeNoNonDuplicateGlobal: " << nodeNoNonDuplicateGlobal << "(removed: " << (nodeNoNonDuplicateGlobal==-1) << ")";

      }
    }
    assert(receiveBufferIndex == nRequestedNodes);
    i++;
  }

  // initialize nNodesLocalWithGhosts_
  nNodesLocalWithGhosts_ = 0;

  // iterate over meshes
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    nNodesLocalWithGhosts_ += subFunctionSpaces_[subMeshIndex]->nNodesLocalWithGhosts();
  }

  nNodesLocalWithGhosts_ -= nNodesSharedLocal_;      // substract double nodes that are included only once
  nNodesLocalWithGhosts_ -= nGhostNodesSharedLocal_; // substract double ghost nodes that are included only once

  LOG(DEBUG) << "meshAndNodeNoLocalToNodeNoNonDuplicateLocal_:   " << meshAndNodeNoLocalToNodeNoNonDuplicateLocal_;
  LOG(DEBUG) << "meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_:  " << meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_;
  LOG(DEBUG) << "isDuplicate_:  " << isDuplicate_;
  LOG(DEBUG) << "nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_:" << nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_;
  LOG(DEBUG) << "nNodesLocalWithGhosts_: " << nNodesLocalWithGhosts_;

}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
createLocalDofOrderings()
{
  MeshPartitionBase::createLocalDofOrderings();

  // fill the dofLocalNo vectors, onlyNodalDofLocalNos_, ghostDofNosGlobalPetsc_ and localToGlobalPetscMappingDofs_
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();

  // initialize onlyNodalDofLocalNos_: vector of local dofs of the nodes, not including derivatives for Hermite
  onlyNodalDofLocalNos_.clear();
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    int nNodesLocalWithoutGhosts = subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts();

    for (dof_no_t nodeNoLocal = 0; nodeNoLocal < nNodesLocalWithoutGhosts; nodeNoLocal++)
    {
      node_no_t nodeNoNonDuplicateLocal = meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal];

      if (nodeNoNonDuplicateLocal != -1)
      {
        dof_no_t dofNoNonDuplicateLocal = nodeNoNonDuplicateLocal*nDofsPerNode;
        onlyNodalDofLocalNos_.push_back(dofNoNonDuplicateLocal);
      }
    }
  }

  // initialize ghostDofNosGlobalPetsc_
  ghostDofNosGlobalPetsc_.clear();
  for (std::vector<PetscInt>::iterator iter = nonDuplicateGhostNodeNosGlobal_.begin(); iter != nonDuplicateGhostNodeNosGlobal_.end(); iter++)
  {
    // loop over dofs on node
    for (int nodalDofNo = 0; nodalDofNo < nDofsPerNode; nodalDofNo++)
    {
      PetscInt nonDuplicateGhostDofNoGlobal = (*iter)*nDofsPerNode + nodalDofNo;
      ghostDofNosGlobalPetsc_.push_back(nonDuplicateGhostDofNoGlobal);
    }
  }
  int nGhostDofs = ghostDofNosGlobalPetsc_.size();

  // create localToGlobalPetscMappingDofs_
  PetscErrorCode ierr;
  Vec temporaryVector;
  ierr = VecCreateGhost(this->mpiCommunicator(), nDofsLocalWithoutGhosts(),
                        nDofsGlobal(), nGhostDofs, ghostDofNosGlobalPetsc_.data(), &temporaryVector); CHKERRV(ierr);

  // retrieve local to global mapping
  ierr = VecGetLocalToGlobalMapping(temporaryVector, &localToGlobalPetscMappingDofs_); CHKERRV(ierr);
  //ierr = VecDestroy(&temporaryVector); CHKERRV(ierr);
  VLOG(1) << "n=" << nDofsLocalWithoutGhosts() << ", N=" << nDofsGlobal() << ", nghost=" << nGhostDofs << " ghosts:" << ghostDofNosGlobalPetsc_;
  VLOG(1) << "Result: " << localToGlobalPetscMappingDofs_;
}

template<int D, typename BasisFunctionType>
std::string MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getString()
{
  std::stringstream result;

  // member variables
  result << "CompositeMesh, nSubMeshes: " << nSubMeshes_
    << ", removedSharedNodes: [";
  int i = 0;
  for (std::vector<std::map<int, std::pair<int,int>>>::iterator iter = removedSharedNodes_.begin(); iter != removedSharedNodes_.end(); iter++, i++)
  {
    for (std::map<int, std::pair<int,int>>::iterator iter2 = iter->begin(); iter2 != iter->end(); iter2++)
    {
      result << "[" << i << "," << iter2->first << "] -> [" << iter2->second.first << "," << iter2->second.second << "] ";
    }
    result << ", ";
  }
  result << "]"
    << ", nElementsLocal: " << nElementsLocal_
    << ", nElementsGlobal: " << nElementsGlobal_
    << ", elementNoGlobalBegin: " << elementNoGlobalBegin_
    << ", nNodesSharedLocal: " << nNodesSharedLocal_
    << ", nGhostNodesSharedLocal: " << nGhostNodesSharedLocal_
    << ", nRemovedNodesNonGhost: [";
  for (int nRemovedNodesNonGhost : nRemovedNodesNonGhost_)
  {
    result << nRemovedNodesNonGhost << ", ";
  }
  result << "]"
    << ", nNonDuplicateNodesWithoutGhosts: [";
  for (int a : nNonDuplicateNodesWithoutGhosts_)
  {
    result << a << ", ";
  }
  result << "]"
    << ", nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts_
    << ", nNodesLocalWithGhosts: " << nNodesLocalWithGhosts_
    << ", nNodesGlobal: " << nNodesGlobal_
    << ", nonDuplicateNodeNoGlobalBegin: " << nonDuplicateNodeNoGlobalBegin_
    << ", meshAndNodeNoLocalToNodeNoNonDuplicateGlobal: [";

  for (std::vector<int> a: meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_)
  {
    result << "[";
    for (int b: a)
      result << b << ",";
    result << "],";
  }
  result << "]"
    << ", meshAndNodeNoLocalToNodeNoNonDuplicateLocal: [";

  for (std::vector<int> a: meshAndNodeNoLocalToNodeNoNonDuplicateLocal_)
  {
    result << "[";
    for (int b: a)
      result << b << ",";
    result << "],";
  }
  result << "]"
    << ", isDuplicate: [";

  for (std::vector<bool> a: isDuplicate_)
  {
    result << "[";
    for (bool b: a)
      result << b << ",";
    result << "],";
  }
  result << "]"
    << ", nodeNoNonDuplicateLocalToMeshAndDuplicateLocal: [";

  for (std::pair<int,node_no_t> a: nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_)
  {
    result << "<" << a.first << "," << a.second << ">,";
  }
  result << "]"
    << ", nonDuplicateGhostNodeNosGlobal: [";
  for (PetscInt a: nonDuplicateGhostNodeNosGlobal_)
    result << a << ",";
  result << "]"
    << ", onlyNodalDofLocalNos: [";
  for (node_no_t a: onlyNodalDofLocalNos_)
    result << a << ",";
  result << "]"
    << ", ghostDofNosGlobalPetsc: [";
  for (PetscInt a: ghostDofNosGlobalPetsc_)
    result << a << ",";
  result << "]";

  // methods
  result << ", nElementsLocal(): " << nElementsLocal()
    << ", nElementsGlobal(): " << nElementsGlobal()
    << ", nDofsLocalWithGhosts(): " << nDofsLocalWithGhosts()
    << ", nDofsLocalWithoutGhosts(): " << nDofsLocalWithoutGhosts()
    << ", nDofsGlobal(): " << nDofsGlobal()
    << ", nNodesLocalWithGhosts(): " << nNodesLocalWithGhosts()
    << ", nNodesLocalWithoutGhosts(): " << nNodesLocalWithoutGhosts()
    << ", nNodesGlobal(): " << nNodesGlobal()
    << ", beginNodeGlobalPetsc(): " << beginNodeGlobalPetsc()
    << ", dofNosLocal(true): [";

  const std::vector<PetscInt> &d0 = dofNosLocal(true);

  for (int i = 0; i < d0.size(); i++)
    result << d0[i] << ",";

  result << "], dofNosLocal(false): [";

  const std::vector<PetscInt> &d1 = dofNosLocal(false);

  for (int i = 0; i < d1.size(); i++)
    result << d1[i] << ",";

  result << "], ghostDofNosGlobalPetsc(): [";

  const std::vector<PetscInt> &d2 = ghostDofNosGlobalPetsc();

  for (int i = 0; i < d2.size(); i++)
    result << d2[i] << ",";

  result << "]";

  // getElementNoGlobalNatural
  result << ", getElementNoGlobalNatural: ";
  for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal(); elementNoLocal++)
  {
    result << elementNoLocal << "->" << getElementNoGlobalNatural(elementNoLocal) << " ";
  }

  // getNodeNoGlobalPetsc
  result << ", getNodeNoGlobalPetsc: ";
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < nNodesLocalWithGhosts(); nodeNoLocal++)
  {
    result << nodeNoLocal << "->" << getNodeNoGlobalPetsc(nodeNoLocal) << " ";
  }

  // getNodeNoGlobalNatural
  result << ", getNodeNoGlobalNatural: ";
  for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal(); elementNoLocal++)
  {
    for (int nodeIndex = 0; nodeIndex < FunctionSpaceType::nNodesPerElement(); nodeIndex++)
    {
      result << elementNoLocal << "," << nodeIndex << "->" << getNodeNoGlobalNatural(elementNoLocal, nodeIndex) << " ";
    }
  }

  // getDofNoGlobalPetsc
  result << ", getDofNoGlobalPetsc: [";
  std::vector<dof_no_t> dofNosLocal(this->nDofsLocalWithGhosts());
  std::iota(dofNosLocal.begin(), dofNosLocal.end(), 0);
  std::vector<PetscInt> dofNosGlobalPetsc;
  getDofNoGlobalPetsc(dofNosLocal, dofNosGlobalPetsc);
  for (PetscInt i : dofNosGlobalPetsc)
    result << i << ",";
  result << "]";

  // getDofNoGlobalPetsc
  result << ", getDofNoGlobalPetsc: ";
  for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithGhosts(); dofNoLocal++)
  {
    result << dofNoLocal << "->" << getDofNoGlobalPetsc(dofNoLocal) << " ";
  }

  // getElementNoGlobalNatural
  result << ", getElementNoLocal: ";
  bool isOnLocalDomain = false;
  for (global_no_t elementNoGlobalPetsc = 0; elementNoGlobalPetsc < nElementsGlobal(); elementNoGlobalPetsc++)
  {
    result << elementNoGlobalPetsc << "->" << getElementNoLocal(elementNoGlobalPetsc, isOnLocalDomain);
    result << "->" << isOnLocalDomain << " ";
  }

  // getNodeNoLocal
  result << ", getNodeNoLocal: ";
  for (global_no_t nodeNoGlobalPetsc = 0; nodeNoGlobalPetsc < nNodesGlobal(); nodeNoGlobalPetsc++)
  {
    result << nodeNoGlobalPetsc << "->" << getNodeNoLocal(nodeNoGlobalPetsc, isOnLocalDomain);
    result << "->" << isOnLocalDomain << " ";
  }

  // getDofNoLocal
  result << ", getDofNoLocal: ";
  for (global_no_t nodeNoGlobalPetsc = 0; nodeNoGlobalPetsc < nNodesGlobal(); nodeNoGlobalPetsc++)
  {
    result << nodeNoGlobalPetsc << "->" << getDofNoLocal(nodeNoGlobalPetsc, isOnLocalDomain);
    result << "->" << isOnLocalDomain << " ";
  }

  // extractLocalNodesWithoutGhosts
  std::vector<global_no_t> globalNodes(this->nNodesGlobal());
  std::iota(globalNodes.begin(), globalNodes.end(), 0);
  extractLocalNodesWithoutGhosts(globalNodes);
  result << ", extractLocalNodesWithoutGhosts: [";
  for (global_no_t i : globalNodes)
    result << i << ",";
  result << "]";

  // extractLocalDofsWithoutGhosts
  std::vector<global_no_t> globalDofs(this->nNodesGlobal());
  std::iota(globalDofs.begin(), globalDofs.end(), 0);
  extractLocalDofsWithoutGhosts(globalDofs);
  result << ", extractLocalDofsWithoutGhosts: [";
  for (global_no_t i : globalDofs)
    result << i << ",";
  result << "]";

  // getSubMeshNoAndElementNoLocal
  result << ", getSubMeshNoAndElementNoLocal: ";
  for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal(); elementNoLocal++)
  {
    int subMeshNo = 0;
    element_no_t elementOnMeshNoLocal = 0;
    getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshNo, elementOnMeshNoLocal);
    result << elementNoLocal << "->" << subMeshNo << "," << elementOnMeshNoLocal << " ";
  }

  // getSubMeshNoAndNodeNoLocal
  result << ", getSubMeshNoAndNodeNoLocal: ";
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < nNodesLocalWithGhosts(); nodeNoLocal++)
  {
    int subMeshNo = 0;
    node_no_t nodeOnMeshNoLocal = 0;
    getSubMeshNoAndNodeNoLocal(nodeNoLocal, subMeshNo, nodeOnMeshNoLocal);
    result << nodeNoLocal << "->" << subMeshNo << "," << nodeOnMeshNoLocal << " ";
  }

  // getSubMeshesWithNodes
  result << ", getSubMeshesWithNodes: ";
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < nNodesLocalWithGhosts(); nodeNoLocal++)
  {
    std::vector<std::pair<int,node_no_t>> subMeshesWithNodes;
    getSubMeshesWithNodes(nodeNoLocal, subMeshesWithNodes);
    result << nodeNoLocal << "->[";
    for (std::pair<int,node_no_t> a : subMeshesWithNodes)
      result << "<" << a.first << "," << a.second << "> ";
    result << "] ";
  }

  return result.str();
}


}  // namespace
