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
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
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

//! get if there are nodes on both borders in the given coordinate direction
//! this is the case if the local partition touches the right/top/back border
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
    // determine if the local partition is at the x+/y+/z+ border of the global domain
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
  ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, dofNosLocal.size(), dofNosLocal.data(), dofNosGlobalPetsc.data()); CHKERRV(ierr);
}

template<typename MeshType,typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNoGlobalPetsc(dof_no_t dofNoLocal) const
{
  PetscInt dofNoGlobal;
  PetscErrorCode ierr;
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
}

template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoLocal(global_no_t nodeNoGlobalPetsc) const
{
  return (node_no_t)(nodeNoGlobalPetsc - nNodesGlobalPetscInPreviousPartitions(ownRankPartitioningIndex_));
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getDofNoLocal(global_no_t dofNoGlobalPetsc) const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  global_no_t nodeNoGlobalPetsc = dofNoGlobalPetsc / nDofsPerNode;
  int nodalDofIndex = dofNoGlobalPetsc % nDofsPerNode;
  return getNodeNoLocal(nodeNoGlobalPetsc) * nDofsPerNode + nodalDofIndex;
}

template<typename MeshType,typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const
{
  std::array<int,MeshType::dim()> coordinatesLocal = getCoordinatesLocal(nodeNoLocal);

  VLOG(2) << "isNonGhost(" << nodeNoLocal << "), coordinatesLocal: " << coordinatesLocal;

  if (MeshType::dim() == 1)
  {
    VLOG(2) << ", nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0) << ")";
  }
  else if (MeshType::dim() == 2)
  {
    VLOG(2) << ", nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1) << ")";
  }
  else if (MeshType::dim() == 3)
  {
    VLOG(2) << ", nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1) << "," << nNodesLocalWithoutGhosts(2) << ")";
  }

  if (nodeNoLocal < nNodesLocalWithoutGhosts())
  {
    return true;
  }

  if (MeshType::dim() == 1)
  {
    if (coordinatesLocal[0] < nNodesLocalWithoutGhosts(0) || hasFullNumberOfNodes(0))    // node is not at x+
    {
      return true;
    }
    else    // node is at x+
    {
      neighbourRankNo = ownRankPartitioningIndex_[0] + 1;
      return false;
    }
  }
  else if (MeshType::dim() == 2)
  {
    VLOG(2) << " isNonGhost(nodeNoLocal=" << nodeNoLocal << "), coordinatesLocal: " << coordinatesLocal
     << ", nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1)
     << ", hasFullNumberOfNodes: " << hasFullNumberOfNodes(0) << "," << hasFullNumberOfNodes(1);

    if (coordinatesLocal[0] < nNodesLocalWithoutGhosts(0) || hasFullNumberOfNodes(0))    // node is not at x+
    {
      if (coordinatesLocal[1] < nNodesLocalWithoutGhosts(1) || hasFullNumberOfNodes(1))    // node is not at y+
      {
        return true;
      }
      else    // node is at y+
      {
        // node is at top (y+)
        neighbourRankNo = (ownRankPartitioningIndex_[1] + 1)*nRanks_[0] + ownRankPartitioningIndex_[0];
        return false;
      }
    }
    else    // node is at x+
    {
      if (coordinatesLocal[1] < nNodesLocalWithoutGhosts(1) || hasFullNumberOfNodes(1))    // node is not at y+
      {
        // node is at right (x+)
        neighbourRankNo = ownRankPartitioningIndex_[1]*nRanks_[0] + ownRankPartitioningIndex_[0] + 1;
        return false;
      }
      else    // node is at y+
      {
        // node is at top right (x+, y+)
        neighbourRankNo = (ownRankPartitioningIndex_[1] + 1)*nRanks_[0] + ownRankPartitioningIndex_[0] + 1;
        return false;
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    if (coordinatesLocal[0] < nNodesLocalWithoutGhosts(0) || hasFullNumberOfNodes(0))    // node is not at x+
    {
      if (coordinatesLocal[1] < nNodesLocalWithoutGhosts(1) || hasFullNumberOfNodes(1))    // node is not at y+
      {
        if (coordinatesLocal[2] < nNodesLocalWithoutGhosts(2) || hasFullNumberOfNodes(2))      // node is not at z+
        {
          return true;
        }
        else        // node is at z+
        {
          // node is at top (z+)
          neighbourRankNo = (ownRankPartitioningIndex_[2] + 1)*nRanks_[0]*nRanks_[1]
            + ownRankPartitioningIndex_[1]*nRanks_[0]
            + ownRankPartitioningIndex_[0];
          VLOG(2) << "node is at top (z+): " << neighbourRankNo;
          return false;
        }
      }
      else   // node is at y+
      {
        if (coordinatesLocal[2] < nNodesLocalWithoutGhosts(2) || hasFullNumberOfNodes(2))    // node is not at z+
        {
          // node is at y+
          neighbourRankNo = ownRankPartitioningIndex_[2]*nRanks_[0]*nRanks_[1]
            + (ownRankPartitioningIndex_[1] + 1)*nRanks_[0]
            + ownRankPartitioningIndex_[0];
          VLOG(2) << "node is at y+: " << neighbourRankNo;
          return false;
        }
        else    // node is at z+
        {
          // node is at y+,z+
          neighbourRankNo = (ownRankPartitioningIndex_[2] + 1)*nRanks_[0]*nRanks_[1]
            + (ownRankPartitioningIndex_[1] + 1)*nRanks_[0]
            + ownRankPartitioningIndex_[0];
          VLOG(2) << "node is at y+,z+: " << neighbourRankNo;
          return false;
        }
      }
    }
    else   // node is at x+
    {
      if (coordinatesLocal[1] < nNodesLocalWithoutGhosts(1) || hasFullNumberOfNodes(1))    // node is not at y+
      {
        if (coordinatesLocal[2] < nNodesLocalWithoutGhosts(2) || hasFullNumberOfNodes(2))      // node is not at z+
        {
          // node is at (x+)
          neighbourRankNo = ownRankPartitioningIndex_[2]*nRanks_[0]*nRanks_[1]
            + ownRankPartitioningIndex_[1]*nRanks_[0]
            + ownRankPartitioningIndex_[0] + 1;
          VLOG(2) << "node is at x+: " << neighbourRankNo;
          return false;
        }
        else      // node is at z+
        {
          // node is at x+,z+
          neighbourRankNo = (ownRankPartitioningIndex_[2] + 1)*nRanks_[0]*nRanks_[1]
            + ownRankPartitioningIndex_[1]*nRanks_[0]
            + ownRankPartitioningIndex_[0] + 1;
          VLOG(2) << "node is at x+,z+: " << neighbourRankNo;
          return false;
        }
      }
      else   // node is at y+
      {
        if (coordinatesLocal[2] < nNodesLocalWithoutGhosts(2) || hasFullNumberOfNodes(2))      // node is not at z+
        {
          // node is at x+,y+
          neighbourRankNo = ownRankPartitioningIndex_[2]*nRanks_[0]*nRanks_[1]
            + (ownRankPartitioningIndex_[1] + 1)*nRanks_[0]
            + (ownRankPartitioningIndex_[0] + 1);
          VLOG(2) << "node is at x+,y+: " << neighbourRankNo;
          return false;
        }
        else      // node is at x+,y+,z+
        {
          // node is at x+,y+,z+
          neighbourRankNo = (ownRankPartitioningIndex_[2] + 1)*nRanks_[0]*nRanks_[1]
            + (ownRankPartitioningIndex_[1] + 1)*nRanks_[0]
            + (ownRankPartitioningIndex_[0] + 1);
          VLOG(2) << "node is at x+,y+,z+: " << neighbourRankNo;
          return false;
        }
      }
    }
  }
  else
  {
    assert(false);
  }
  return false;  // this is only needed for cray compiler, it will not be reached
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
neighbourRank(Mesh::face_t face)
{
  /*
  face0Minus = 0, face0Plus,
  face1Minus, face1Plus,
  face2Minus, face2Plus
  * */
  if (face == Mesh::face_t::face0Minus)
  {
    // if this subdomain is at the left end of the global domain
    if (ownRankPartitioningIndex_[0] == 0)
    {
      return -1;
    }
    else
    {
      return this->ownRankNo()-1;
    }
  }
  else if (face == Mesh::face_t::face0Plus)
  {
    // if this subdomain is at the right end of the global domain
    if (ownRankPartitioningIndex_[0] == nRanks_[0]-1)
    {
      return -1;
    }
    else
    {
      return this->ownRankNo()+1;
    }
  }
  else if (face == Mesh::face_t::face1Minus)
  {
    assert(MeshType::dim() >= 2);

    // if this subdomain is at the front end of the global domain
    if (ownRankPartitioningIndex_[1] == 0)
    {
      return -1;
    }
    else
    {
      int neighbourRankNo = (ownRankPartitioningIndex_[1]-1)*nRanks_[0] + ownRankPartitioningIndex_[0];  // 2D case
      if (MeshType::dim() == 3)
      {
        neighbourRankNo += ownRankPartitioningIndex_[2]*nRanks_[0]*nRanks_[1];
      }
      return neighbourRankNo;
    }
  }
  else if (face == Mesh::face_t::face1Plus)
  {
    assert(MeshType::dim() >= 2);

    // if this subdomain is at the back end of the global domain
    if (ownRankPartitioningIndex_[1] == nRanks_[1]-1)
    {
      return -1;
    }
    else
    {
      int neighbourRankNo = (ownRankPartitioningIndex_[1] + 1)*nRanks_[0]
        + ownRankPartitioningIndex_[0];  // 2D case
      if (MeshType::dim() == 3)
      {
        neighbourRankNo += ownRankPartitioningIndex_[2]*nRanks_[0]*nRanks_[1];
      }
      return neighbourRankNo;
    }
  }
  else if (face == Mesh::face_t::face2Minus)
  {
    assert(MeshType::dim() == 3);

    // if this subdomain is at the bottom end of the global domain
    if (ownRankPartitioningIndex_[2] == 0)
    {
      return -1;
    }
    else
    {
      return (ownRankPartitioningIndex_[2] - 1)*nRanks_[0]*nRanks_[1]
        + ownRankPartitioningIndex_[1]*nRanks_[0]
        + ownRankPartitioningIndex_[0];
    }
  }
  else if (face == Mesh::face_t::face2Plus)
  {
    assert(MeshType::dim() == 3);

    // if this subdomain is at the top end of the global domain
    if (ownRankPartitioningIndex_[2] == nRanks_[2]-1)
    {
      return -1;
    }
    else
    {
      return (ownRankPartitioningIndex_[2] + 1)*nRanks_[0]*nRanks_[1]
        + ownRankPartitioningIndex_[1]*nRanks_[0]
        + ownRankPartitioningIndex_[0];
    }
  }
  return -1;  // does not happen (but intel compiler does not recognize it)
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,MeshType::dim()> &nBoundaryElements, std::vector<dof_no_t> &dofNos)
{
  neighbourRankNo = neighbourRank(face);

  // if there is no neighbouring rank, do not fill data structures
  if (neighbourRankNo == -1)
  {
    return;
  }

  /*
  face0Minus = 0, face0Plus,
  face1Minus, face1Plus,
  face2Minus, face2Plus
   */

  std::array<element_no_t,MeshType::dim()> boundaryElementIndexStart;

  if (face == Mesh::face_t::face0Minus || face == Mesh::face_t::face0Plus)
  {
    if (face == Mesh::face_t::face0Minus)
    {
      boundaryElementIndexStart[0] = 0;
      nBoundaryElements[0] = 1;
    }
    else
    {
      boundaryElementIndexStart[0] = nElementsLocal_[0]-1;
      nBoundaryElements[0] = 1;
    }

    if (MeshType::dim() >= 2)
    {
      boundaryElementIndexStart[1] = 0;
      nBoundaryElements[1] = nElementsLocal_[1];
    }
    if (MeshType::dim() == 3)
    {
      boundaryElementIndexStart[2] = 0;
      nBoundaryElements[2] = nElementsLocal_[2];
    }
  }
  else if (face == Mesh::face_t::face1Minus || face == Mesh::face_t::face1Plus)
  {
    assert(MeshType::dim() >= 2);

    if (face == Mesh::face_t::face1Minus)
    {
      boundaryElementIndexStart[1] = 0;
      nBoundaryElements[1] = 1;
    }
    else
    {
      boundaryElementIndexStart[1] = nElementsLocal_[1]-1;
      nBoundaryElements[1] = 1;
    }

    boundaryElementIndexStart[0] = 0;
    nBoundaryElements[0] = nElementsLocal_[0];

    if (MeshType::dim() == 3)
    {
      boundaryElementIndexStart[2] = 0;
      nBoundaryElements[2] = nElementsLocal_[2];
    }
  }
  else if (face == Mesh::face_t::face2Minus || face == Mesh::face_t::face2Plus)
  {
    assert(MeshType::dim() == 3);

    if (face == Mesh::face_t::face2Minus)
    {
      boundaryElementIndexStart[2] = 0;
      nBoundaryElements[2] = 1;
    }
    else
    {
      boundaryElementIndexStart[2] = nElementsLocal_[2]-1;
      nBoundaryElements[2] = 1;
    }

    boundaryElementIndexStart[0] = 0;
    nBoundaryElements[0] = nElementsLocal_[0];

    boundaryElementIndexStart[1] = 0;
    nBoundaryElements[1] = nElementsLocal_[1];
  }

  const int averageNDofsPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
  const int nDofsPerNode = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode();

  int nBoundaryDofsTotal = 1;
  switch (MeshType::dim())
  {
  case 1:
    nBoundaryDofsTotal = (nBoundaryElements[0]*averageNDofsPerElement1D + nDofsPerNode);
    break;
  case 2:
    nBoundaryDofsTotal = (nBoundaryElements[0]*averageNDofsPerElement1D + nDofsPerNode) * (nBoundaryElements[1]*averageNDofsPerElement1D + nDofsPerNode);
    break;
  case 3:
    nBoundaryDofsTotal = (nBoundaryElements[0]*averageNDofsPerElement1D + nDofsPerNode) * (nBoundaryElements[1]*averageNDofsPerElement1D + nDofsPerNode) * (nBoundaryElements[2]*averageNDofsPerElement1D + nDofsPerNode);
    break;
  }

  LOG(DEBUG) << "getBoundaryElements: nBoundaryElements: " << nBoundaryElements[0] << "x" << nBoundaryElements[1] << "x" << nBoundaryElements[2]
    << ", nBoundaryDofsTotal: " << nBoundaryDofsTotal
    << " nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1) << "," << nNodesLocalWithoutGhosts(2) << ")"
    << ", nNodesLocalWithGhosts: (" << nNodesLocalWithGhosts(0) << "," << nNodesLocalWithGhosts(1) << "," << nNodesLocalWithGhosts(2) << ")";

  // determine dofs of all nodes adjacent to the boundary elements
  dofNos.resize(nBoundaryDofsTotal);


  std::array<dof_no_t,MeshType::dim()> nDofs;
  std::array<dof_no_t,MeshType::dim()> boundaryDofIndexStart;
  for (int dimensionIndex = 0; dimensionIndex < MeshType::dim(); dimensionIndex++)
  {
    nDofs[dimensionIndex] = nBoundaryElements[dimensionIndex]*averageNDofsPerElement1D + nDofsPerNode;
    boundaryDofIndexStart[dimensionIndex] = boundaryElementIndexStart[dimensionIndex]*averageNDofsPerElement1D;
  }

  LOG(DEBUG) << "nDofs: " << nDofs << ", boundaryElementIndexStart: " << boundaryElementIndexStart << ", boundaryDofIndexStart: " << boundaryDofIndexStart;

  if (MeshType::dim() == 1)
  {
    int dofIndex = 0;
    for (int i = boundaryDofIndexStart[0]; i < boundaryDofIndexStart[0]+nDofs[0]; i++, dofIndex++)
    {
      dofNos[dofIndex] = dofNosLocalNaturalOrdering_[i];
    }
  }
  else if (MeshType::dim() == 2)
  {
    int nDofsRow = nElementsLocal_[0]*averageNDofsPerElement1D + averageNDofsPerElement1D;
    int dofIndex = 0;
    for (int j = boundaryDofIndexStart[1]; j < boundaryDofIndexStart[1]+nDofs[1]; j++)
    {
      for (int i = boundaryDofIndexStart[0]; i < boundaryDofIndexStart[0]+nDofs[0]; i++, dofIndex++)
      {
        dofNos[dofIndex] = dofNosLocalNaturalOrdering_[j*nDofsRow + i];
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    int nDofsRow = nElementsLocal_[0]*averageNDofsPerElement1D + averageNDofsPerElement1D;
    int nDofsPlane = (nElementsLocal_[1]*averageNDofsPerElement1D + averageNDofsPerElement1D) * nDofsRow;
    //LOG(DEBUG) << "nDofsRow: " << nDofsRow << ", nDofsPlane: " << nDofsPlane;
    int dofIndex = 0;
    for (int k = boundaryDofIndexStart[2]; k < boundaryDofIndexStart[2]+nDofs[2]; k++)
    {
      for (int j = boundaryDofIndexStart[1]; j < boundaryDofIndexStart[1]+nDofs[1]; j++)
      {
        for (int i = boundaryDofIndexStart[0]; i < boundaryDofIndexStart[0]+nDofs[0]; i++, dofIndex++)
        {
          dofNos[dofIndex] = dofNosLocalNaturalOrdering_[k*nDofsPlane + j*nDofsRow + i];
          //LOG(DEBUG) << "   (i,j,k) = (" << i << "," << j << "," << k << "), dofIndex=" << dofIndex << "<" << nBoundaryDofsTotal << ", index " << k*nDofsPlane + j*nDofsRow + i << "<" << dofNosLocalNaturalOrdering_.size()
          //  << ", dofNo: " << dofNos[dofIndex];

        }
      }
    }
  }
}

template<typename MeshType,typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
ownRankPartitioningIndex(int coordinateDirection)
{
  return ownRankPartitioningIndex_[coordinateDirection];
}

}  // namespace
