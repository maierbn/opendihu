#include "basis_on_mesh/04_basis_on_mesh_numbers_structured.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

namespace BasisOnMesh
{

// These numberings run over the whole locally stored information, i.e. there is no distinguishing between ghost and interior nodes here.
 
// element-local dofIndex to local dofNo for 1D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement:
  // 1         2            2
  VLOG(3) << "getDofNo<1D>(elementNo=" << elementNo << ", dofIndex=" << dofIndex << ") = " << BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
  
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 1D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local dofIndex to local dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;
  
  return getNodeNo(elementNo, nodeIndex)*nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 2D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local dofIndex to local dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;
  
  return getNodeNo(elementNo, nodeIndex)*nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 3D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local nodeIndex to local nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
node_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0 1
  // averageNDofsPerElement:
  // 1         2            2
  // nDofsPerBasis:
  // 2         3            4
  // nDofsPerNode:
  // 1         1            2
  // nNodesPerElement:
  // 2         3            2
  
  VLOG(3) << "getNodeNo<1D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " 
    << BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
  
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
}

// element-local nodeIndex to local nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
node_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4

  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  
  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionLocal_;
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  
  // the number of non-ghost nodes in different rows
  node_no_t nodesPerRow = averageNNodesPerElement1D * nElements[0] + (this->meshPartition()->hasFullNumberOfNodes(0)? 1: 0);
  element_no_t elementX = element_no_t(elementNo % nElements[0]);
  element_no_t elementY = element_no_t(elementNo / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  // check for ghost nodes which have nos after the normal contiguous numbering scheme
  if (!this->meshPartition()->hasFullNumberOfNodes(1) && elementY == nElements[1]-1 && localY == nNodesPerElement1D-1)  
  {
    // node is a ghost node on the top border   
    // if there are ghost nodes on the right border
    if (this->meshPartition()->hasFullNumberOfNodes(0))
    {
      VLOG(3) << "getNodeNo<2D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " 
        << this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * nElements[1] + averageNNodesPerElement1D * elementX + localX;
      
      return this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * nElements[1]
        + averageNNodesPerElement1D * elementX + localX;
    }
    else 
    {
      VLOG(3) << "getNodeNo<2D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " 
        << this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * elementX + localX;
      
      return this->meshPartition()->nNodesLocalWithoutGhosts()
        + averageNNodesPerElement1D * elementX + localX;
    }
  }
  
  if (!this->meshPartition()->hasFullNumberOfNodes(0) && elementX == nElements[0]-1 && localX == nNodesPerElement1D-1)
  {
    VLOG(3) << "getNodeNo<2D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " 
      << this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * elementY + localY;
    
    // node is a ghost node on the right border
    return this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * elementY + localY;
  }
  
  VLOG(3) << "getNodeNo<2D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " 
    << nodesPerRow * (elementY * averageNNodesPerElement1D + localY) + averageNNodesPerElement1D * elementX + localX;
    
  // compute local node no for non-ghost node
  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex to local nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
node_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given

  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionLocal_;
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  node_no_t nodesPerRow0 = (averageNNodesPerElement1D * nElements[0] + (this->meshPartition()->hasFullNumberOfNodes(0)? 1: 0));
  node_no_t nodesPerPlane = (averageNNodesPerElement1D * nElements[1] + (this->meshPartition()->hasFullNumberOfNodes(1)? 1: 0)) * nodesPerRow0;

  element_no_t elementZ = element_no_t(elementNo / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  element_no_t elementX = elementNo % nElements[0];
  dof_no_t localZ = dof_no_t(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  dof_no_t localY = dof_no_t((nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;

  // check for ghost nodes which have nos after the normal contiguous numbering scheme
  // handle ghosts on z+ border
  if (!this->meshPartition()->hasFullNumberOfNodes(2) && elementZ == nElements[2]-1 && localZ == nNodesPerElement1D-1)
  {
    // node is a ghost node on the z+ border, i.e. lies in a x-y plane on top of the non-ghost dofs
    node_no_t nodeNo = this->meshPartition()->nNodesLocalWithoutGhosts();
    if (this->meshPartition()->hasFullNumberOfNodes(0))
    {
      if (!this->meshPartition()->hasFullNumberOfNodes(1))
      {
        nodeNo += nodesPerRow0 * (elementZ * averageNNodesPerElement1D + localZ);
      }
    }
    else 
    {
      if (this->meshPartition()->hasFullNumberOfNodes(1))
      {
        node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1] + 1;
        nodeNo += nodesPerRow1 * (elementZ * averageNNodesPerElement1D + localZ);
      }
      else 
      {
        node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1];
        nodeNo += (nodesPerRow0 + nodesPerRow1 + 1) * (elementZ * averageNNodesPerElement1D + localZ);
      }
    }
    
    nodeNo += (averageNNodesPerElement1D * nElements[0] + 1) * (elementY * averageNNodesPerElement1D + localY);
    nodeNo += averageNNodesPerElement1D * elementX + localX;
    
    VLOG(3) << "getNodeNo<3D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " << nodeNo;
    return nodeNo;
  }
  
  // handle ghosts on y+ border
  if (!this->meshPartition()->hasFullNumberOfNodes(1) && elementY == nElements[1]-1 && localY == nNodesPerElement1D-1)
  {
    // node is a ghost node on the y+ border, i.e. lies in a x-z plane
    node_no_t nodeNo = this->meshPartition()->nNodesLocalWithoutGhosts();
    if (this->meshPartition()->hasFullNumberOfNodes(0))
    {
      nodeNo += nodesPerRow0 * (elementZ * averageNNodesPerElement1D + localZ);
    }
    else 
    {
      node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1];
      nodeNo += (nodesPerRow0 + nodesPerRow1 + 1) * (elementZ * averageNNodesPerElement1D + localZ);
      nodeNo += nodesPerRow1;
    }
    
    nodeNo += averageNNodesPerElement1D * elementX + localX;
    
    VLOG(3) << "getNodeNo<3D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " << nodeNo;
    return nodeNo;
  }
  
  // handle ghosts on x+ border
  if (!this->meshPartition()->hasFullNumberOfNodes(0) && elementX == nElements[0]-1 && localX == nNodesPerElement1D-1)
  {
    // node is a ghost node on the x+ border, i.e. lies in a y-z plane
    node_no_t nodeNo = this->meshPartition()->nNodesLocalWithoutGhosts();
    if (this->meshPartition()->hasFullNumberOfNodes(1))
    {
      node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1] + 1;
      nodeNo += nodesPerRow1 * (elementZ * averageNNodesPerElement1D + localZ);
    }
    else 
    {
      node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1];
      nodeNo += (nodesPerRow0 + nodesPerRow1 + 1) * (elementZ * averageNNodesPerElement1D + localZ);
    }
    nodeNo += averageNNodesPerElement1D * elementY + localY;
    
    VLOG(3) << "getNodeNo<3D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " << nodeNo;
    return nodeNo;
  }
  
  VLOG(3) << "getNodeNo<3D>(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = " 
    << nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ) + nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY) + averageNNodesPerElement1D * elementX + localX;
    
  // compute local node no for non-ghost node
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ)
    + nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex of global element to global nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
global_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNoGlobal(global_no_t elementNoGlobal, int nodeIndex) const
{
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNoGlobal + nodeIndex;
}

// element-local nodeIndex of global element to global nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
global_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeNoGlobal(global_no_t elementNoGlobal, int nodeIndex) const
{
  const std::array<global_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionGlobal_;
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();

  // the number of non-ghost nodes in different rows
  global_no_t nodesPerRow = averageNNodesPerElement1D * nElements[0] + 1;
  element_no_t elementX = element_no_t(elementNoGlobal % nElements[0]);
  element_no_t elementY = element_no_t(elementNoGlobal / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex of global element to global nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
global_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNoGlobal(global_no_t elementNoGlobal, int nodeIndex) const
{
  const std::array<global_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionGlobal_;
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  global_no_t nodesPerRow0 = averageNNodesPerElement1D * nElements[0] + 1;
  global_no_t nodesPerPlane = (averageNNodesPerElement1D * nElements[1] + 1) * nodesPerRow0;

  element_no_t elementZ = element_no_t(elementNoGlobal / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t((elementNoGlobal % (nElements[0] * nElements[1])) / nElements[0]);
  element_no_t elementX = elementNoGlobal % nElements[0];
  dof_no_t localZ = dof_no_t(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  dof_no_t localY = dof_no_t((nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;

  // compute local node no for non-ghost node
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ)
    + nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

};  // namespace
