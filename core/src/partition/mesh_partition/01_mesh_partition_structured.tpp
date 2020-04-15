#include "partition/mesh_partition/01_mesh_partition_structured.h"

#include <cstdlib>
#include "utility/vector_operators.h"
#include "function_space/00_function_space_base_dim.h"

namespace Partition
{

template<typename MeshType,typename BasisFunctionType>
MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
MeshPartition(std::array<global_no_t,MeshType::dim()> nElementsGlobal, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), nElementsGlobal_(nElementsGlobal), isDegenerate_(false), hasFullNumberOfNodes_({false}), nDofsLocalWithoutGhosts_(-1)
{
  VLOG(1) << "create MeshPartition where only the global size is known, " 
    << "nElementsGlobal: " << nElementsGlobal_ << ", rankSubset: " << *rankSubset << ", mesh dimension: " << MeshType::dim();
  
  if (MeshType::dim() == 1 && nElementsGlobal_[0] == 0)
  {
    initialize1NodeMesh();
  }
  else if (!rankSubset->ownRankIsContained())
  {
    isDegenerate_ = true;
    initializeDegenerateMesh();
  }
  else
  {
    // determine partitioning of elements
    this->createDmElements();

    // initialize the cached value of nDofsLocalWithoutGhosts
    this->setNDofsLocalWithoutGhosts();

    // initialize dof vectors
    this->createLocalDofOrderings();
  }
  
  LOG(DEBUG) << "nElementsLocal_: " << nElementsLocal_ << ", nElementsGlobal_: " << nElementsGlobal_
    << ", hasFullNumberOfNodes_: " << hasFullNumberOfNodes_;
  LOG(DEBUG) << "nRanks: " << nRanks_ << ", localSizesOnPartitions_: " << localSizesOnPartitions_
    << ", beginElementGlobal_: " << beginElementGlobal_;
  
  LOG(DEBUG) << *this;
}

template<typename MeshType,typename BasisFunctionType>
MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
MeshPartition(std::array<node_no_t,MeshType::dim()> nElementsLocal, std::array<global_no_t,MeshType::dim()> nElementsGlobal, 
              std::array<global_no_t,MeshType::dim()> beginElementGlobal, 
              std::array<int,MeshType::dim()> nRanks, std::shared_ptr<RankSubset> rankSubset) :
  MeshPartitionBase(rankSubset), beginElementGlobal_(beginElementGlobal), nElementsLocal_(nElementsLocal), nElementsGlobal_(nElementsGlobal), 
  nRanks_(nRanks), isDegenerate_(false), hasFullNumberOfNodes_({false}), nDofsLocalWithoutGhosts_(-1)
{
  // partitioning is already prescribed as every rank knows its own local size
 
  LOG(DEBUG) << "create MeshPartition where every rank already knows its own local size. "
    << "nElementsLocal: " << nElementsLocal << ", nElementsGlobal: " << nElementsGlobal 
    << ", beginElementGlobal: " << beginElementGlobal << ", nRanks: " << nRanks << ", rankSubset: " << *rankSubset;
    
  if (MeshType::dim() == 1 && nElementsGlobal_[0] == 0)
  {
    initialize1NodeMesh();
  }
  else if (!rankSubset->ownRankIsContained())
  {
    isDegenerate_ = true;
    initializeDegenerateMesh();
  }
  else
  {
    initializeHasFullNumberOfNodes();

    // determine localSizesOnRanks
    std::array<std::vector<element_no_t>,MeshType::dim()> localSizesOnRanks;    // [dimensionIndex][rankNo]
    for (int i = 0; i < MeshType::dim(); i++)
    {
      localSizesOnRanks[i].resize(rankSubset->size());
    }

    for (int i = 0; i < MeshType::dim(); i++)
    {
      MPIUtility::handleReturnValue(MPI_Allgather(&nElementsLocal_[i], 1, MPIU_INT,
        localSizesOnRanks[i].data(), 1, MPIU_INT, rankSubset->mpiCommunicator()));
    }
    LOG(DEBUG) << "determined localSizesOnRanks: " << localSizesOnRanks;
    LOG(DEBUG) << "MeshType::dim(): " << MeshType::dim() << ", nRanks: " << nRanks_;
    
    // create localSizesOnPartitions_ from localSizesOnRanks
    for (int dimensionIndex = 0; dimensionIndex < MeshType::dim(); dimensionIndex++)
    {
      VLOG(1) << "dimensionIndex: " << dimensionIndex << ", resize to " << nRanks_[dimensionIndex];
      //assert (nRanks_[dimensionIndex] != 0);
      localSizesOnPartitions_[dimensionIndex].resize(nRanks_[dimensionIndex]);

      // loop over the first rank of the respective portion
      int rankStride = 1;
      if (dimensionIndex == 0)
      {
        rankStride = 1;
      }
      else if (dimensionIndex == 1)
      {
        rankStride = nRanks_[0];
      }
      else if (dimensionIndex == 2)
      {
        rankStride = nRanks_[0]*nRanks_[1];
      }

      VLOG(1) << "   rankStride: " << rankStride;
      int partitionIndex = 0;
      for (int rankNo = 0; partitionIndex < nRanks_[dimensionIndex]; rankNo += rankStride, partitionIndex++)
      {
        VLOG(1) << "   rankNo: " << rankNo;
        VLOG(1) << "   localSizesOnRanks: " << localSizesOnRanks[dimensionIndex][rankNo];
        localSizesOnPartitions_[dimensionIndex][partitionIndex] = localSizesOnRanks[dimensionIndex][rankNo];
        VLOG(1) << "    saved to partitionIndex " << partitionIndex;
        if (partitionIndex >= nRanks_[dimensionIndex]-1)
          break;
      }
    }

    LOG(DEBUG) << "determined localSizesOnPartitions_: " << localSizesOnPartitions_;

    setOwnRankPartitioningIndex();

    // initialize the cached value of nDofsLocalWithoutGhosts
    this->setNDofsLocalWithoutGhosts();

    this->createLocalDofOrderings();
  }

  LOG(DEBUG) << *this;

  for (int i = 0; i < MeshType::dim(); i++)
  {
    LOG(DEBUG) << "  beginNodeGlobalNatural(" << i << "): " << beginNodeGlobalNatural(i);
  }
}

//! get the partition index in a given coordinate direction from the rankNo
template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
convertRankNoToPartitionIndex(int coordinateDirection, int rankNo)
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  //direction 0 is fastest, then 1, then 2
  if (coordinateDirection == 0)
  {
    return rankNo % nRanks_[0];
  }
  else if (coordinateDirection == 1)
  {
    return (rankNo % (nRanks_[0]*nRanks_[1])) / nRanks_[0];
  }
  else if (coordinateDirection == 2)
  {
    return rankNo / (nRanks_[0]*nRanks_[1]);
  }

  return 0;   // will not be reached
}
  
template<typename MeshType,typename BasisFunctionType>
template<typename T>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents) const
{
  LOG(DEBUG) << "extractLocalNodesWithoutGhosts, nComponents=" << nComponents << ", input: " << vector;

  std::vector<T> result(nNodesLocalWithoutGhosts()*nComponents);
  global_no_t resultIndex = 0;
  
  if (MeshType::dim() == 1)
  {
    assert(vector.size() >= (beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0)) * nComponents);
    for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
    {
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        result[resultIndex++] = vector[i*nComponents + componentNo];
      }
    }
  }
  else if (MeshType::dim() == 2)
  {
    assert(vector.size() >= ((beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) - 1)*nNodesGlobal(0)
      + beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0)) * nComponents);
    for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1); j++)
    {
      for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
      {
        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          result[resultIndex++] = vector[(j*nNodesGlobal(0) + i)*nComponents + componentNo];
        }
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    assert(vector.size() >= ((beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2) - 1)*nNodesGlobal(1)*nNodesGlobal(0)
      + (beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) - 1)*nNodesGlobal(0)
      + beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0)) * nComponents);
    for (global_no_t k = beginNodeGlobalNatural(2); k < beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2); k++)
    {
      for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1); j++)
      {
        for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
        {
          for (int componentNo = 0; componentNo < nComponents; componentNo++)
          {
            result[resultIndex++] = vector[(k*nNodesGlobal(1)*nNodesGlobal(0) + j*nNodesGlobal(0) + i)*nComponents + componentNo];
          }
        }
      }
    }
  }
  LOG(DEBUG) << "                   output: " << result;
  
  // store values
  vector.assign(result.begin(), result.end());
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
extractLocalDofsWithoutGhosts(std::vector<double> &vector) const
{
  this->template extractLocalDofsWithoutGhosts<double>(vector);
}
  
template<typename MeshType,typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
extractLocalDofsWithoutGhosts(std::vector<T> &vector) const
{
  dof_no_t nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  std::vector<T> result(nDofsLocalWithoutGhosts());
  global_no_t resultIndex = 0;
  
  if (MeshType::dim() == 1)
  {
    if (vector.size() < (beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0))*nDofsPerNode)
    {
      LOG(DEBUG) << "vector.size: " << vector.size() << ", beginNodeGlobalNatural(0): " << beginNodeGlobalNatural(0)
        << ", nNodesLocalWithoutGhosts(0): " << nNodesLocalWithoutGhosts(0)
        << ", nDofsPerNode: " << nDofsPerNode << *this;
    }
    assert(vector.size() >= (beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0))*nDofsPerNode);
    for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
    {
      for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
      {
        result[resultIndex++] = vector[i*nDofsPerNode + dofOnNodeIndex];
      }
    }
  }
  else if (MeshType::dim() == 2)
  {
    assert(vector.size() >= ((beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) - 1)*nNodesGlobal(0) + (beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0)))*nDofsPerNode);
    for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1); j++)
    {
      for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
      {
        for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
        {
          result[resultIndex++] = vector[(j*nNodesGlobal(0) + i)*nDofsPerNode + dofOnNodeIndex];
        }
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    if (vector.size() < ((beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2)-1)*nNodesGlobal(1)*nNodesGlobal(0)
      + (beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1)-1)*nNodesGlobal(0) + beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0))*nDofsPerNode)
    {
      LOG(FATAL) << "vector.size(): " << vector.size() << ", expected: " << ((beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2)-1)*nNodesGlobal(1)*nNodesGlobal(0)
        + (beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1)-1)*nNodesGlobal(0) + beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0))*nDofsPerNode
        << ", nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1) << "," << nNodesLocalWithoutGhosts(2)
        << ", nNodesGlobal: " << nNodesGlobal(0) << "," << nNodesGlobal(1) << "," << nNodesGlobal(2) << ", nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts();        
    }

    assert(vector.size() >= ((beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2)-1)*nNodesGlobal(1)*nNodesGlobal(0)
      + (beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1)-1)*nNodesGlobal(0) + beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0))*nDofsPerNode);
    for (global_no_t k = beginNodeGlobalNatural(2); k < beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2); k++)
    {
      for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1); j++)
      {
        for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
        {
          for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
          {
            result[resultIndex++] = vector[(k*nNodesGlobal(1)*nNodesGlobal(0) + j*nNodesGlobal(0) + i)*nDofsPerNode + dofOnNodeIndex];
          }
        }
      }
    }
  }
  
  assert(resultIndex == nDofsLocalWithoutGhosts());

  // store values
  vector.assign(result.begin(), result.end());
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
refine(std::array<int,MeshType::dim()> refinementFactor)
{
  LOG(DEBUG) << "refine mesh partition by refinementFactors " << refinementFactor;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    beginElementGlobal_[i] *= refinementFactor[i];   ///< global element no.s of the lower left front corner of the domain
    nElementsLocal_[i] *= refinementFactor[i];     ///< local size, i.e. number of nodes in the coordinate directions of the local portion
    nElementsGlobal_[i] *= refinementFactor[i];    ///< global number of elements in the coodinate directions

    ///< the sizes of different partitions in each coordinate direction, i.e. localSizesOnPartitions_[0] is (width partition #0, width partition #1, ...)
    for (int partitionIndex = 0; partitionIndex < localSizesOnPartitions_[i].size(); partitionIndex++)
    {
      localSizesOnPartitions_[i][partitionIndex] *= refinementFactor[i];
    }
  }

  setNDofsLocalWithoutGhosts();

  //! fill the dofLocalNo vectors, onlyNodalDofLocalNos_, ghostDofNosGlobalPetsc_ and localToGlobalPetscMappingDofs_
  createLocalDofOrderings();

  dofNosLocalNaturalOrdering_.clear();

  LOG(DEBUG) << "mesh partition after refinement: " << *this;
}



}  // namespace
