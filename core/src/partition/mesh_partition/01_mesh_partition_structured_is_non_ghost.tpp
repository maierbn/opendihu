#include "partition/mesh_partition/01_mesh_partition.h"

#include <cstdlib>
#include "utility/vector_operators.h"
#include "function_space/00_function_space_base_dim.h"

namespace Partition
{

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

}  // namespace
