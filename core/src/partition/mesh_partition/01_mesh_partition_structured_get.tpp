#include "partition/mesh_partition/01_mesh_partition.h"

#include <cstdlib>
#include "utility/vector_operators.h"
#include "function_space/00_function_space_base_dim.h"

namespace Partition
{


//! get the local to global mapping for the current partition
template<typename MeshType,typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
localToGlobalMappingDofs()
{
  return localToGlobalPetscMappingDofs_;
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nRanks(int coordinateDirection) const
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());
  return nRanks_[coordinateDirection];
}

//! number of entries in the current partition
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nElementsLocal() const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  VLOG(1) << "determine nElementsLocal";
  element_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    VLOG(1) << "     * " << nElementsLocal_[i];
    result *= nElementsLocal_[i];
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nElementsGlobal() const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  global_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    result *= nElementsGlobal_[i];
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nDofsLocalWithGhosts() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();

  return nNodesLocalWithGhosts() * nDofsPerNode;
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nDofsLocalWithoutGhosts() const
{
  return this->nDofsLocalWithoutGhosts_;
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nDofsGlobal() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();

  return nNodesGlobal() * nDofsPerNode;
}

//! number of nodes in the local partition
template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesLocalWithGhosts() const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  element_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    result *= nNodesLocalWithGhosts(i);
  }
  return result;
}

//! number of nodes in the local partition
template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesLocalWithoutGhosts() const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  element_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    result *= nNodesLocalWithoutGhosts(i);
  }
  return result;
}

//! number of nodes in the local partition
template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesLocalWithGhosts(int coordinateDirection, int partitionIndex) const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  const int nNodesPer1DElement = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();

  if (partitionIndex == -1)
  {
    //VLOG(2) << "nNodesLocalWithGhosts(coordinateDirection=" << coordinateDirection << ", partitionIndex=" << partitionIndex
    //  << "), nElementsLocal_: " << nElementsLocal_ << ", nNodesPer1DElement: " << nNodesPer1DElement
    //  << ", result: " << this->nElementsLocal_[coordinateDirection] * nNodesPer1DElement + 1;

    // if partitionIndex was given as -1, it means return the the value for the own partition
    return this->nElementsLocal_[coordinateDirection] * nNodesPer1DElement + 1;
  }
  else
  {
    //VLOG(2) << "nNodesLocalWithGhosts(coordinateDirection=" << coordinateDirection << ", partitionIndex=" << partitionIndex
    //  << "), localSizesOnPartitions: " << localSizesOnPartitions_ << ", nNodesPer1DElement: " << nNodesPer1DElement
    //  << ", result: " << this->localSizesOnPartitions_[coordinateDirection][partitionIndex] * nNodesPer1DElement + 1;

    // get the value for the given partition with index partitionIndex
    return this->localSizesOnPartitions_[coordinateDirection][partitionIndex] * nNodesPer1DElement + 1;
  }
}

//! number of nodes in the local partition
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nElementsLocal(int coordinateDirection) const
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  return this->nElementsLocal_[coordinateDirection];
}

//! number of nodes in total
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nElementsGlobal(int coordinateDirection) const
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  return this->nElementsGlobal_[coordinateDirection];
}

//! number of nodes in the local partition
template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesLocalWithoutGhosts(int coordinateDirection, int partitionIndex) const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  const int nNodesPer1DElement = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();

  if (partitionIndex == -1)
  {
    // if partitionIndex was given as -1, it means return the the value for the own partition
    return this->nElementsLocal_[coordinateDirection] * nNodesPer1DElement
      + (this->hasFullNumberOfNodes(coordinateDirection)? 1 : 0);
  }
  else
  {
    // get the value for the given partition with index partitionIndex
    //VLOG(1) << "nNodesLocalWithoutGhosts(coordinateDirection=" << coordinateDirection << ", partitionIndex=" << partitionIndex << "), localSizesOnPartitions: "
    //  << this->localSizesOnPartitions_ << ": " << this->localSizesOnPartitions_[coordinateDirection][partitionIndex] << ", has full: " << this->hasFullNumberOfNodes(coordinateDirection, partitionIndex);
    return this->localSizesOnPartitions_[coordinateDirection][partitionIndex] * nNodesPer1DElement
      + (this->hasFullNumberOfNodes(coordinateDirection, partitionIndex)? 1 : 0);
  }
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
beginElementGlobal(int coordinateDirection) const
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  return beginElementGlobal_[coordinateDirection];
}

//! number of nodes in total
template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesGlobal() const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  global_no_t result = 1;
  for (int i = 0; i < MeshType::dim(); i++)
  {
    result *= nNodesGlobal(i);
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
beginNodeGlobalNatural(int coordinateDirection, int partitionIndex) const
{
  //for degenerate mesh
  if (isDegenerate_)
    return 0;

  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  const int nNodesPer1DElement = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();

  if (partitionIndex == -1)
  {
    //VLOG(2) << "beginNodeGlobalNatural(coordinateDirection=" << coordinateDirection << ", partitionIndex=" << partitionIndex
    //  << "), beginElementGlobal: " << beginElementGlobal_ << ", result: " << beginElementGlobal_[coordinateDirection] * nNodesPer1DElement;

    // if partitionIndex was given as -1, it means return the global natural no of the beginNode for the current partition
    return beginElementGlobal_[coordinateDirection] * nNodesPer1DElement;
  }
  else
  {

    // compute the global natural no of the beginNode for partition with partitionIndex in coordinateDirection
    global_no_t nodeNoGlobalNatural = 0;
    // loop over all partitions in the given coordinateDirection that are before the partition with partitionIndex
    for (int i = 0; i < partitionIndex; i++)
    {
      // sum up the number of nodes on these previous partitions
      nodeNoGlobalNatural += (localSizesOnPartitions_[coordinateDirection][i]) * nNodesPer1DElement;
    }

    //VLOG(2) << "beginNodeGlobalNatural(coordinateDirection=" << coordinateDirection << ", partitionIndex=" << partitionIndex
    //  << "), localSizesOnPartitions_: " << localSizesOnPartitions_ << ", result: " << nodeNoGlobalNatural;

    return nodeNoGlobalNatural;
  }
}

//! number of nodes in total
template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesGlobal(int coordinateDirection) const
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  return nElementsGlobal_[coordinateDirection] * FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement() + 1;
}

//! get if there are nodes on both boundarys in the given coordinate direction
//! this is the case if the local partition touches the right/top/back boundary
template<typename MeshType,typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
hasFullNumberOfNodes(int coordinateDirection, int partitionIndex) const
{
  assert(0 <= coordinateDirection);
  assert(coordinateDirection < MeshType::dim());

  if (partitionIndex == -1)
  {
    // if partitionIndex was given as -1, it means return the information for the own partition
    return hasFullNumberOfNodes_[coordinateDirection];
  }
  else
  {
    // determine if the local partition is at the x+/y+/z+ boundary of the global domain
    return (partitionIndex == nRanks_[coordinateDirection]-1);
  }
}

//! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
//! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
template<typename MeshType,typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
dofNosLocal(bool onlyNodalValues) const
{
  if (onlyNodalValues)
  {
    return onlyNodalDofLocalNos_;
  }
  else
  {
    return this->dofNosLocal_;
  }
}

template<typename MeshType,typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
ghostDofNosGlobalPetsc() const
{
  return ghostDofNosGlobalPetsc_;
}

template<typename MeshType,typename BasisFunctionType>
const std::vector<dof_no_t> &MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
dofNosLocalNaturalOrdering() const
{
  return dofNosLocalNaturalOrdering_;
}

//! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const
{
  dof_no_t nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  dofNosGlobalNatural.resize(nDofsLocalWithoutGhosts());
  global_no_t resultIndex = 0;

  VLOG(1) << "getDofNosGlobalNatural, nDofsLocal: " << nDofsLocalWithoutGhosts();

  if (MeshType::dim() == 1)
  {
    for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
    {
      for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
      {
        dofNosGlobalNatural[resultIndex++] = i*nDofsPerNode + dofOnNodeIndex;
      }
    }
  }
  else if (MeshType::dim() == 2)
  {
    for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1); j++)
    {
      for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
      {
        for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
        {
          dofNosGlobalNatural[resultIndex++] = (j*nNodesGlobal(0) + i)*nDofsPerNode + dofOnNodeIndex;
        }
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    for (global_no_t k = beginNodeGlobalNatural(2); k < beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2); k++)
    {
      for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1); j++)
      {
        for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i++)
        {
          for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
          {
            dofNosGlobalNatural[resultIndex++] = (k*nNodesGlobal(1)*nNodesGlobal(0) + j*nNodesGlobal(0) + i)*nDofsPerNode + dofOnNodeIndex;
          }
        }
      }
    }
  }
  else
  {
    assert(false);
  }
}

template<typename MeshType,typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getPartitioningIndex(std::array<global_no_t,MeshType::dim()> nodeNoGlobalNatural) const
{
  std::array<int,MeshType::dim()> partitioningIndex;

  const int nNodesPer1DElement = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();

  // determine the x-index in the partition grid, i.e. how many partitions (other ranks) there are in x direction left of the current ranks' partition
  int partitionX = 0;
  global_no_t xGlobalNatural = 0;
  while (xGlobalNatural <= nodeNoGlobalNatural[0] && xGlobalNatural < nNodesGlobal(0)-1)
  {
    VLOG(3) << "   x GlobalNatural=" << xGlobalNatural << ", partitionX=" << partitionX << ", nodeNoGlobalNatural[0]=" << nodeNoGlobalNatural[0];
    assert(localSizesOnPartitions_[0].size() > partitionX);
    xGlobalNatural += localSizesOnPartitions_[0][partitionX++]*nNodesPer1DElement;
    VLOG(3) << "   x GlobalNatural=" << xGlobalNatural << ", partitionX=" << partitionX << ", nodeNoGlobalNatural[0]=" << nodeNoGlobalNatural[0];
  }
  partitionX--;

  // store partitionX
  partitioningIndex[0] = partitionX;

  if (MeshType::dim() >= 2)
  {
    // determine the y-index in the partition grid, i.e. how many partitions (other ranks) there are in y direction before the current ranks' partition
    int partitionY = 0;
    global_no_t yGlobalNatural = 0;
    while (yGlobalNatural <= nodeNoGlobalNatural[1] && yGlobalNatural < nNodesGlobal(1)-1)
    {
      assert(localSizesOnPartitions_[1].size() > partitionY);
      yGlobalNatural += localSizesOnPartitions_[1][partitionY++]*nNodesPer1DElement;
      VLOG(3) << "   y GlobalNatural=" << yGlobalNatural << ", partitionY=" << partitionY << ", nodeNoGlobalNatural[1]=" << nodeNoGlobalNatural[1];
    }
    partitionY--;

    // store partitionY
    partitioningIndex[1] = partitionY;
  }

  if (MeshType::dim() == 3)
  {
    // determine the z-index in the partition grid, i.e. how many partitions (other ranks) there are in z direction before the current ranks' partition
    int partitionZ = 0;
    global_no_t zGlobalNatural = 0;
    while (zGlobalNatural <= nodeNoGlobalNatural[2] && zGlobalNatural < nNodesGlobal(2)-1)
    {
      assert(localSizesOnPartitions_[2].size() > partitionZ);
      zGlobalNatural += localSizesOnPartitions_[2][partitionZ++]*nNodesPer1DElement;
      VLOG(3) << "   z GlobalNatural=" << zGlobalNatural << ", partitionZ=" << partitionZ << ", nodeNoGlobalNatural[2]=" << nodeNoGlobalNatural[2];
    }
    partitionZ--;

    // store partitionZ
    partitioningIndex[2] = partitionZ;
  }

  return partitioningIndex;
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
beginNodeGlobalPetsc() const
{
  return nNodesGlobalPetscInPreviousPartitions(ownRankPartitioningIndex_);
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
nNodesGlobalPetscInPreviousPartitions(std::array<int,MeshType::dim()> partitionIndex) const
{
  if (MeshType::dim() == 1)
  {
    return beginNodeGlobalNatural(0,partitionIndex[0]);
  }
  else if (MeshType::dim() == 2)
  {
    //  _ _   globalNatural nos in the partitions before current partition
    // |_|_L_    <- (2)
    // |_|_|_|   <- (1)
    // |_|_|_|

    VLOG(3) << "beginNodeGlobalNatural: " << beginNodeGlobalNatural(0,partitionIndex[0]) << " " << beginNodeGlobalNatural(1,partitionIndex[1]);
    VLOG(3) << "nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts(0,partitionIndex[0]) << " " << nNodesLocalWithoutGhosts(1,partitionIndex[1]);

    // compute globalNatural no from globalPetsc no (i,j)
    return beginNodeGlobalNatural(1,partitionIndex[1])*nNodesGlobal(0)  // (1)
      + beginNodeGlobalNatural(0,partitionIndex[0])*nNodesLocalWithoutGhosts(1,partitionIndex[1]);   // (2)
  }
  else if (MeshType::dim() == 3)
  {
    VLOG(3) << "beginNodeGlobalNatural: " << beginNodeGlobalNatural(0,partitionIndex[0]) << " " << beginNodeGlobalNatural(1,partitionIndex[1]) << " " << beginNodeGlobalNatural(2,partitionIndex[2]);
    VLOG(3) << "nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts(0,partitionIndex[0]) << " " << nNodesLocalWithoutGhosts(1,partitionIndex[1]) << " " << nNodesLocalWithoutGhosts(2,partitionIndex[2]);

    return beginNodeGlobalNatural(2,partitionIndex[2])*nNodesGlobal(1)*nNodesGlobal(0)
      + nNodesLocalWithoutGhosts(2,partitionIndex[2])*beginNodeGlobalNatural(1,partitionIndex[1])*nNodesGlobal(0)  // (1)
      + nNodesLocalWithoutGhosts(2,partitionIndex[2])*nNodesLocalWithoutGhosts(1,partitionIndex[1])*beginNodeGlobalNatural(0,partitionIndex[0]);   //(2)

  }
  else
  {
    assert(false);
  }
  return 0;  // this never happens, but the intel compiler does not recognize it (gcc does)
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const
{
  std::array<global_no_t,MeshType::dim()> coordinatesGlobal = getCoordinatesGlobal(nodeNoLocal);
  return getNodeNoGlobalPetsc(coordinatesGlobal);
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getRankOfNodeNoGlobalNatural(global_no_t nodeNoGlobalNatural) const
{
  // determine global coordinates of node
  std::array<global_no_t,MeshType::dim()> coordinatesGlobal = getCoordinatesGlobalOfNodeNoGlobalNatural(nodeNoGlobalNatural);

  // determine partitionIndex
  std::array<int,MeshType::dim()> partitionIndex = getPartitioningIndex(coordinatesGlobal);

  int rankNo = partitionIndex[0];
  if (MeshType::dim() == 2)
  {
    rankNo = partitionIndex[1]*nRanks_[0] + partitionIndex[0];
  }
  else if (MeshType::dim() == 3)
  {
    rankNo = partitionIndex[2]*nRanks_[0]*nRanks_[1] + partitionIndex[1]*nRanks_[0] + partitionIndex[0];
  }

  return rankNo;
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getRankOfDofNoGlobalNatural(global_no_t dofNoGlobalNatural) const
{
  global_no_t nodeNoGlobalNatural = dofNoGlobalNatural / FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  return getRankOfNodeNoGlobalNatural(nodeNoGlobalNatural);
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoGlobalPetsc(std::array<global_no_t,MeshType::dim()> coordinatesGlobal) const
{
  std::array<int,MeshType::dim()> partitionIndex = getPartitioningIndex(coordinatesGlobal);

  if (MeshType::dim() == 1)
  {
    global_no_t nodeNoGlobalPetsc = coordinatesGlobal[0];
    return nodeNoGlobalPetsc;
  }
  else if (MeshType::dim() == 2)
  {
    // compute globalNatural no from globalPetsc no (i,j)
    global_no_t nodeNoGlobalPetsc = nNodesGlobalPetscInPreviousPartitions(partitionIndex)
      + (coordinatesGlobal[1]-beginNodeGlobalNatural(1,partitionIndex[1])) * nNodesLocalWithoutGhosts(0,partitionIndex[0])
      + (coordinatesGlobal[0]-beginNodeGlobalNatural(0,partitionIndex[0]));   // in-partition local no
    return nodeNoGlobalPetsc;
  }
  else if (MeshType::dim() == 3)
  {
    // compute globalNatural no from globalPetsc no (i,j,k)
    global_no_t nodeNoGlobalPetsc = nNodesGlobalPetscInPreviousPartitions(partitionIndex)
      + (coordinatesGlobal[2]-beginNodeGlobalNatural(2,partitionIndex[2]))*nNodesLocalWithoutGhosts(1,partitionIndex[1])*nNodesLocalWithoutGhosts(0,partitionIndex[0])
      + (coordinatesGlobal[1]-beginNodeGlobalNatural(1,partitionIndex[1]))*nNodesLocalWithoutGhosts(0,partitionIndex[0])
      + (coordinatesGlobal[0]-beginNodeGlobalNatural(0,partitionIndex[0]));   // in-partition local no

    VLOG(1) << "coordinates global: " << coordinatesGlobal << ", partitionIndex: " << partitionIndex
      << ", nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0,partitionIndex[0]) << "," << nNodesLocalWithoutGhosts(1,partitionIndex[1]) << "," << nNodesLocalWithoutGhosts(2,partitionIndex[2]) << ")"
      << ", beginNodeGlobalNatural: (" << beginNodeGlobalNatural(0,partitionIndex[0]) << "," << beginNodeGlobalNatural(1,partitionIndex[1]) << "," << beginNodeGlobalNatural(2,partitionIndex[2]) << "), "
      << "nodeNoGlobalPetsc: " << nodeNoGlobalPetsc
      << " = " << nNodesGlobalPetscInPreviousPartitions(partitionIndex) << ", ";

    return nodeNoGlobalPetsc;
  }
  else
  {
    assert(false);
  }
  return 0;  // this never happens, but the intel compiler does not recognize it (gcc does)
}

//! get the node no in global petsc ordering from a local node no
template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const
{
  dofNosGlobalPetsc.resize(dofNosLocal.size());

  // transfer the local indices to global indices
  PetscErrorCode ierr;
  if (localToGlobalPetscMappingDofs_)
    ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, dofNosLocal.size(), dofNosLocal.data(), dofNosGlobalPetsc.data()); CHKERRV(ierr);
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNoGlobalPetsc(dof_no_t dofNoLocal) const
{
  assert(dofNoLocal >= 0);
  assert(dofNoLocal < nDofsLocalWithGhosts());

  PetscInt dofNoGlobal;
  PetscErrorCode ierr;
  if (localToGlobalPetscMappingDofs_)
    ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, 1, &dofNoLocal, &dofNoGlobal); CHKERRQ(ierr);
  return (global_no_t)dofNoGlobal;
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getElementNoGlobalNatural(element_no_t elementNoLocal) const
{
  if (MeshType::dim() == 1)
  {
    return beginElementGlobal_[0] + elementNoLocal;
  }
  else if (MeshType::dim() == 2)
  {
    element_no_t elementY = element_no_t(elementNoLocal / nElementsLocal(0));
    element_no_t elementX = elementNoLocal % nElementsLocal(0);

    return (beginElementGlobal_[1] + elementY) * nElementsGlobal_[0]
      + beginElementGlobal_[0] + elementX;
  }
  else if (MeshType::dim() == 3)
  {
    element_no_t elementZ = element_no_t(elementNoLocal / (nElementsLocal(0) * nElementsLocal(1)));
    element_no_t elementY = element_no_t((elementNoLocal % (nElementsLocal(0) * nElementsLocal(1))) / nElementsLocal(0));
    element_no_t elementX = elementNoLocal % nElementsLocal(0);

    return (beginElementGlobal_[2] + elementZ) * nElementsGlobal_[0] * nElementsGlobal_[1]
      + (beginElementGlobal_[1] + elementY) * nElementsGlobal_[0]
      + beginElementGlobal_[0] + elementX;
  }
  else
  {
    assert(false);
  }
  return 0;
}


  //! get the local element no. from coordinates
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getElementNoLocal(std::array<int,MeshType::dim()> elementCoordinates) const
{
  if (MeshType::dim() == 1)
  {
    return elementCoordinates[0];
  }
  else if (MeshType::dim() == 2)
  {
    return elementCoordinates[1]*nElementsLocal_[0] + elementCoordinates[0];
  }
  else if (MeshType::dim() == 3)
  {
    return elementCoordinates[2]*nElementsLocal_[0]*nElementsLocal_[1] + elementCoordinates[1]*nElementsLocal_[0] + elementCoordinates[0];
  }
  else
  {
    assert(false);
  }
  return 0;
}

//! get the local element no. from the global no., set isOnLocalDomain to true if the node with global coordinates is in the local domain
template<typename MeshType,typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const
{
  // get global coordinates of the lower front left node of the element
  std::array<global_no_t,MeshType::dim()> coordinatesGlobal;
  coordinatesGlobal[0] = elementNoGlobalPetsc % nElementsGlobal(0);

  if (MeshType::dim() >= 2)
  {
    coordinatesGlobal[1] = (elementNoGlobalPetsc % (nElementsGlobal(0)*nElementsGlobal(1))) / nElementsGlobal(0);
  }

  if (MeshType::dim() >= 3)
  {
    coordinatesGlobal[2] = elementNoGlobalPetsc / (nElementsGlobal(0)*nElementsGlobal(1));
  }

  // get local coordinates
  std::array<int,MeshType::dim()> coordinatesLocal = getCoordinatesLocal(coordinatesGlobal, isOnLocalDomain);

  if (isOnLocalDomain)
  {
    // get local element no
    element_no_t elementNoLocal = getElementNoLocal(coordinatesLocal);
    return elementNoLocal;
  }

  return 0;
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoGlobalNatural(std::array<global_no_t,MeshType::dim()> coordinatesGlobal) const
{
  if (MeshType::dim() == 1)
  {
    return coordinatesGlobal[0];
  }
  else if (MeshType::dim() == 2)
  {
    return coordinatesGlobal[1]*nNodesGlobal(0) + coordinatesGlobal[0];
  }
  else if (MeshType::dim() == 3)
  {
    return coordinatesGlobal[2]*nNodesGlobal(0)*nNodesGlobal(1) + coordinatesGlobal[1]*nNodesGlobal(0) + coordinatesGlobal[0];
  }
  else
  {
    assert(false);
  }
  return 0; // never reached, but pgi compiler does not recognize it
}

template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const
{
  node_no_t nodeNoLocal = (node_no_t)(nodeNoGlobalPetsc - nNodesGlobalPetscInPreviousPartitions(ownRankPartitioningIndex_));
  isLocal = (nodeNoLocal >= 0) && (nodeNoLocal < nDofsLocalWithoutGhosts());
  return nodeNoLocal;
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  global_no_t nodeNoGlobalPetsc = dofNoGlobalPetsc / nDofsPerNode;
  int nodalDofIndex = dofNoGlobalPetsc % nDofsPerNode;
  return getNodeNoLocal(nodeNoGlobalPetsc, isLocal) * nDofsPerNode + nodalDofIndex;
}

template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoLocalFromGlobalNatural(global_no_t nodeNoGlobalNatural, bool &isOnLocalDomain) const
{
  std::array<global_no_t,MeshType::dim()> coordinatesGlobal = getCoordinatesGlobalOfNodeNoGlobalNatural(nodeNoGlobalNatural);
  return getNodeNoLocal(coordinatesGlobal, isOnLocalDomain);
}

//! get the local node no for its global coordinates
template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, bool &isOnLocalDomain) const
{
  std::array<int,MeshType::dim()> coordinatesLocal = getCoordinatesLocal(coordinatesGlobal, isOnLocalDomain);

  if (isOnLocalDomain)
  {
    node_no_t nodeNoLocal = this->getNodeNoLocal(coordinatesLocal);
    return nodeNoLocal;
  }
  return -1;
}

//! get the local dof no for the global coordinates of the node
template<typename MeshType,typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNoLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, int nodalDofIndex, bool &isOnLocalDomain) const
{
  node_no_t nodeNoLocal = getNodeNoLocal(coordinatesGlobal, isOnLocalDomain);
  if (isOnLocalDomain)
    return nodeNoLocal * FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode() + nodalDofIndex;

  return -1;
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
ownRankPartitioningIndex(int coordinateDirection)
{
  return ownRankPartitioningIndex_[coordinateDirection];
}

//! check if the element is at an outer corner in the x-y plane of the domain, then set isAtCorner to true and return which corner (one of edge0Minus1Minus, edge0Plus1Minus, edge0Minus1Plus,  edge0Plus1Plus), only for 3D meshes!
template<typename MeshType,typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
elementIsAtCorner(element_no_t elementNoLocal, Mesh::face_or_edge_t &edge)
{
  element_no_t elementX = elementNoLocal % nElementsLocal(0);
  element_no_t elementY = element_no_t((elementNoLocal % (nElementsLocal(0) * nElementsLocal(1))) / nElementsLocal(0));

  if (ownRankPartitioningIndex_[0] == 0 && elementX == 0)   // left
  {
    if (ownRankPartitioningIndex_[1] == 0 && elementY == 0) // front
    {
      // left front
      edge = Mesh::face_or_edge_t::edge0Minus1Minus;
      return true;
    }
    else if (ownRankPartitioningIndex_[1] == nRanks_[1]-1 && elementY == nElementsGlobal_[1]-1)  // back
    {
      // left back
      edge = Mesh::face_or_edge_t::edge0Minus1Plus;
      return true;
    }
  }
  else if (ownRankPartitioningIndex_[0] == nRanks_[0]-1 && elementX == nElementsGlobal_[0]-1)   // right
  {
    if (ownRankPartitioningIndex_[1] == 0 && elementY == 0)
    {
      // right front
      edge = Mesh::face_or_edge_t::edge0Plus1Minus;
      return true;
    }
    else if (ownRankPartitioningIndex_[1] == nRanks_[1]-1 && elementY == nElementsGlobal_[1]-1)   // back
    {
      // right back
      edge = Mesh::face_or_edge_t::edge0Plus1Plus;
      return true;
    }
  }
  return false;
}


}  // namespace
