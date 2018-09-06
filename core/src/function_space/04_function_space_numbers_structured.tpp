#include "function_space/04_function_space_numbers_structured.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

namespace FunctionSpace
{

// These numberings run over the whole locally stored information, i.e. there is no distinguishing between ghost and interior nodes here.
 
// element-local dofIndex to local dofNo for 1D
template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement:
  // 1         2            2
  VLOG(3) << "getDofNo<1D>(elementNo=" << elementNo << ", dofIndex=" << dofIndex << ") = " << FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
  
  return FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 1D
template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local dofIndex to local dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;
  
  return getNodeNo(elementNo, nodeIndex)*nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 2D
template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local dofIndex to local dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;
  
  return getNodeNo(elementNo, nodeIndex)*nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 3D
template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local nodeIndex to local nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
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
    << FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
  
  return FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
}

// element-local nodeIndex to local nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
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
  int averageNNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  
  // the number of non-ghost nodes in different rows
  node_no_t nodesPerRow = averageNNodesPerElement1D * nElements[0] + (this->meshPartition()->hasFullNumberOfNodes(0)? 1: 0);
  element_no_t elementX = element_no_t(elementNo % nElements[0]);
  element_no_t elementY = element_no_t(elementNo / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  VLOG(3) << "elementY=" << elementY << ", elementX=" << elementX;
  VLOG(3) << "localX=" << localX << ", localY=" << localY;

  // check for ghost nodes which have nos after the normal contiguous numbering scheme
  if (!this->meshPartition()->hasFullNumberOfNodes(1) && elementY == nElements[1]-1 && localY == nNodesPerElement1D-1)  
  {
    // node is a ghost node on the top border   
    // if there are ghost nodes on the right border
    if (!this->meshPartition()->hasFullNumberOfNodes(0))
    {
      VLOG(3) << "getNodeNo<2D>b(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = "
        << this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * nElements[1] + averageNNodesPerElement1D * elementX + localX;
      
      return this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * nElements[1]
        + averageNNodesPerElement1D * elementX + localX;
    }
    else 
    {
      VLOG(3) << "getNodeNo<2D>c(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = "
        << this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * elementX + localX;
      
      return this->meshPartition()->nNodesLocalWithoutGhosts()
        + averageNNodesPerElement1D * elementX + localX;
    }
  }
  
  if (!this->meshPartition()->hasFullNumberOfNodes(0) && elementX == nElements[0]-1 && localX == nNodesPerElement1D-1)
  {
    VLOG(3) << "getNodeNo<2D>d(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = "
      << this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * elementY + localY;
    
    // node is a ghost node on the right border
    return this->meshPartition()->nNodesLocalWithoutGhosts() + averageNNodesPerElement1D * elementY + localY;
  }
  
  VLOG(3) << "getNodeNo<2D>a(elementNo=" << elementNo << ", nodeIndex=" << nodeIndex << ") = "
    << nodesPerRow * (elementY * averageNNodesPerElement1D + localY) + averageNNodesPerElement1D * elementX + localX;
    
  // compute local node no for non-ghost node
  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex to local nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given

  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionLocal_;
  int averageNNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
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
global_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
{
  return FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNoGlobalNatural + nodeIndex;
}

// element-local nodeIndex of global element to global nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
global_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
{
  const std::array<global_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionGlobal_;
  int averageNNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();

  // the number of non-ghost nodes in different rows
  global_no_t nodesPerRow = averageNNodesPerElement1D * nElements[0] + 1;
  element_no_t elementX = element_no_t(elementNoGlobalNatural % nElements[0]);
  element_no_t elementY = element_no_t(elementNoGlobalNatural / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex of global element to global nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
global_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
{
  const std::array<global_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirectionGlobal_;
  int averageNNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  global_no_t nodesPerRow0 = averageNNodesPerElement1D * nElements[0] + 1;
  global_no_t nodesPerPlane = (averageNNodesPerElement1D * nElements[1] + 1) * nodesPerRow0;

  element_no_t elementZ = element_no_t(elementNoGlobalNatural / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t((elementNoGlobalNatural % (nElements[0] * nElements[1])) / nElements[0]);
  element_no_t elementX = elementNoGlobalNatural % nElements[0];
  node_no_t localZ = node_no_t(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  node_no_t localY = node_no_t((nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  node_no_t localX = nodeIndex % nNodesPerElement1D;

  // compute local node no for non-ghost node
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ)
    + nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// local node no of neighbour node, may be a ghost node, for 1D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>>::
getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const
{
  assert (nodeNoLocal < this->nNodesLocalWithoutGhosts());

  if (direction == Mesh::face_t::face0Minus)  // left node
  {
    // if there is no neighbouring node
    if (nodeNoLocal == 0)
    {
      return -1;
    }

    // get previous node no
    return nodeNoLocal - 1;
  }
  else if (direction == Mesh::face_t::face0Plus) // right node
  {
    // if there is no neighbouring node
    if (nodeNoLocal == this->meshPartition()->nNodesLocalWithGhosts(0)-1)
    {
      return -1;
    }

    // get next node no, this is correct for ghost as well as non-ghost nodes
    return nodeNoLocal + 1;
  }
  else
  {
    assert(false);
  }
}

// local node no of neighbour node, may be a ghost node, for 2D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>>::
getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const
{
  assert (nodeNoLocal < this->nNodesLocalWithoutGhosts());

  dof_no_t localX = nodeNoLocal % this->meshPartition()->nNodesLocalWithoutGhosts(0);
  dof_no_t localY = dof_no_t(nodeNoLocal / this->meshPartition()->nNodesLocalWithoutGhosts(0));


  if (direction == Mesh::face_t::face0Minus)  // left node
  {
    // if there is no neighbouring node
    if (localX == 0)
    {
      return -1;
    }
    // get previous node no
    return nodeNoLocal - 1;
  }
  else if (direction == Mesh::face_t::face0Plus) // right node
  {
    // if the node is on the right row, there is no neighbouring node
    if (localX == this->meshPartition()->nNodesLocalWithGhosts(0)-1)
    {
      return -1;
    }

    // check if right node is a ghost node
    if (localX == this->meshPartition()->nNodesLocalWithoutGhosts(0)-1)
    {
      // there is a row of ghost nodes on the right
      return this->meshPartition()->nNodesLocalWithoutGhosts() + localY;
    }

    // return next node no
    return nodeNoLocal + 1;
  }
  else if (direction == Mesh::face_t::face1Minus)  // bottom node
  {
    // if there is no neighbouring node
    if (localY == 0)
    {
      return -1;
    }

    // get node below
    return nodeNoLocal - this->meshPartition()->nNodesLocalWithoutGhosts(0);
  }
  else if (direction == Mesh::face_t::face1Plus)  // top node
  {
    // if the node is at the top row, there is no neighbouring node
    if (localY == this->meshPartition()->nNodesLocalWithGhosts(1)-1)
    {
      return -1;
    }

    if (localY == this->meshPartition()->nNodesLocalWithoutGhosts(1)-1)
    {
      // there is one row of ghost nodes above the current
      if (this->meshPartition()->hasFullNumberOfNodes(0))
      {
        // there are only top ghost nodes
        return this->meshPartition()->nNodesLocalWithoutGhosts() + localX;
      }
      else
      {
        // there are right and top ghost nodes
        return this->meshPartition()->nNodesLocalWithoutGhosts() + this->meshPartition()->nNodesLocalWithGhosts(1) - 1 + localX;
      }
    }

    // get node above
    return nodeNoLocal + this->meshPartition()->nNodesLocalWithoutGhosts(0);
  }
  else
  {
    assert(false);
  }
}

// local node no of neighbour node, may be a ghost node, for 3D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>>::
getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const
{
  assert (nodeNoLocal < this->nNodesLocalWithoutGhosts());

  node_no_t localZ = node_no_t(nodeNoLocal / (this->meshPartition()->nNodesLocalWithoutGhosts(0) * this->meshPartition()->nNodesLocalWithoutGhosts(1)));
  node_no_t localY = node_no_t((nodeNoLocal % (this->meshPartition()->nNodesLocalWithoutGhosts(0) * this->meshPartition()->nNodesLocalWithoutGhosts(1))) / this->meshPartition()->nNodesLocalWithoutGhosts(0));
  node_no_t localX = nodeNoLocal % this->meshPartition()->nNodesLocalWithoutGhosts(0);

  if (direction == Mesh::face_t::face0Minus)  // left node
  {
    // if there is no neighbouring node
    if (localX == 0)
    {
      return -1;
    }

    // get previous node no
    return nodeNoLocal - 1;
  }
  else if (direction == Mesh::face_t::face0Plus) // right node
  {
    // if the node is on the right row, there is no neighbouring node
    if (localX == this->meshPartition()->nNodesLocalWithGhosts(0)-1)
    {
      return -1;
    }

    // if right node is a ghost node, this implies !this->meshPartition()->hasFullNumberOfNodes(0)
    if (localX == this->meshPartition()->nNodesLocalWithoutGhosts(0)-1)
    {
      // there is a y-z plane of ghost nodes on the right
      node_no_t neighbourNodeNo = this->meshPartition()->nNodesLocalWithoutGhosts();

      if (this->meshPartition()->hasFullNumberOfNodes(1))
      {
        neighbourNodeNo += (this->meshPartition()->nNodesLocalWithGhosts(1)) * localZ;
      }
      else
      {
        neighbourNodeNo += (this->meshPartition()->nNodesLocalWithGhosts(0) + this->meshPartition()->nNodesLocalWithGhosts(1) - 1) * localZ;
      }
      return neighbourNodeNo + localY;
    }

    // return next node no
    return nodeNoLocal + 1;
  }
  else if (direction == Mesh::face_t::face1Minus)  // y- node
  {
    // if there is no neighbouring node
    if (localY == 0)
    {
      return -1;
    }

    // get node below
    return nodeNoLocal - this->meshPartition()->nNodesLocalWithoutGhosts(0);
  }
  else if (direction == Mesh::face_t::face1Plus)  // y+ node
  {
    LOG(DEBUG) << "direction Y+, localX=" << localX << ", localY=" << localY << ", localZ=" << localZ << ", pf" << this->meshPartition()->nNodesLocalWithoutGhosts(1)-1;

    // if the node is at the top row, there is no neighbouring node
    if (localY == this->meshPartition()->nNodesLocalWithGhosts(1)-1)
    {
      return -1;
    }

    // if y+ node is a ghost node, this implies !this->meshPartition()->hasFullNumberOfNodes(1)
    if (localY == this->meshPartition()->nNodesLocalWithoutGhosts(1)-1)
    {
      // there is one row of ghost nodes behind the current
      node_no_t neighbourNodeNo = this->meshPartition()->nNodesLocalWithoutGhosts();

      if (this->meshPartition()->hasFullNumberOfNodes(0))
      {
        // there are only ghosts at y+
        neighbourNodeNo += (this->meshPartition()->nNodesLocalWithoutGhosts(0)) * localZ;
      }
      else
      {
        // there are ghosts at y+ and x+
        neighbourNodeNo += (this->meshPartition()->nNodesLocalWithGhosts(0) + this->meshPartition()->nNodesLocalWithGhosts(1) - 1) * localZ;
        neighbourNodeNo += this->meshPartition()->nNodesLocalWithoutGhosts(1);
      }
      return neighbourNodeNo + localX;
    }

    // get node behind
    return nodeNoLocal + this->meshPartition()->nNodesLocalWithoutGhosts(0);
  }
  else if (direction == Mesh::face_t::face2Minus)  // z- node
  {
    // if there is no node below
    if (localZ == 0)
    {
      return -1;
    }

    // get node below
    return nodeNoLocal - this->meshPartition()->nNodesLocalWithoutGhosts(0)*this->meshPartition()->nNodesLocalWithoutGhosts(1);
  }
  else if (direction == Mesh::face_t::face2Plus)  // z+ node
  {
    // if the node is at the top row, there is no neighbouring node
    if (localZ == this->meshPartition()->nNodesLocalWithGhosts(2)-1)
    {
      return -1;
    }

    // if z+ node is a ghost node, this implies !this->meshPartition()->hasFullNumberOfNodes(2)
    if (localZ == this->meshPartition()->nNodesLocalWithoutGhosts(2)-1)
    {
      // there is one row of ghost nodes above the current
      node_no_t neighbourNodeNo = this->meshPartition()->nNodesLocalWithoutGhosts();

      LOG(DEBUG) << "direction Z+, localX=" << localX << ", localY=" << localY << ", localZ=" << localZ << ", starting at " << neighbourNodeNo
        << ", nNodesLocalWithoutGhosts: " << this->meshPartition()->nNodesLocalWithoutGhosts(0) << "," << this->meshPartition()->nNodesLocalWithoutGhosts(1) << "," << this->meshPartition()->nNodesLocalWithoutGhosts(2)
        << ", hasFullNumberOfNodes: " << this->meshPartition()->hasFullNumberOfNodes(0) << "," << this->meshPartition()->hasFullNumberOfNodes(1) << "," << this->meshPartition()->hasFullNumberOfNodes(2);

      if (this->meshPartition()->hasFullNumberOfNodes(0))
      {
        if (!this->meshPartition()->hasFullNumberOfNodes(1))
        {
          // there are ghosts at y+ and z+
          LOG(DEBUG) << "a";
          neighbourNodeNo += this->meshPartition()->nNodesLocalWithoutGhosts(0) * this->meshPartition()->nNodesLocalWithoutGhosts(2);
        }
      }
      else
      {
        if (this->meshPartition()->hasFullNumberOfNodes(1))
        {
          // there are ghosts at x+ and z+
          LOG(DEBUG) << "B";
          neighbourNodeNo += this->meshPartition()->nNodesLocalWithoutGhosts(1) * this->meshPartition()->nNodesLocalWithoutGhosts(2);
        }
        else
        {
          // there are ghosts at x+, y+ and z+
          LOG(DEBUG) << "C";
          neighbourNodeNo += (this->meshPartition()->nNodesLocalWithGhosts(0) + this->meshPartition()->nNodesLocalWithGhosts(1) - 1) * this->meshPartition()->nNodesLocalWithoutGhosts(2);
        }
      }

      if (this->meshPartition()->hasFullNumberOfNodes(0))
      {
        LOG(DEBUG) << "D";
        // there are no ghosts at x+
        return neighbourNodeNo + this->meshPartition()->nNodesLocalWithoutGhosts(0)*localY + localX;
      }
      else
      {
        LOG(DEBUG) << "E";
        // there are ghosts at x+
        return neighbourNodeNo + this->meshPartition()->nNodesLocalWithGhosts(0)*localY + localX;
      }
    }

    // get node above
    return nodeNoLocal + this->meshPartition()->nNodesLocalWithoutGhosts(0)*this->meshPartition()->nNodesLocalWithoutGhosts(1);
  }
  else
  {
    assert(false);
  }
}

// get local node no from the coordinates (x), may be a ghost node, for 1D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>>::
getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const
{
  return coordinateLocal[0];
}

// get local node no from the coordinates (x,y), may be a ghost node, for 2D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>>::
getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const
{
  dof_no_t localX = coordinateLocal[0];
  dof_no_t localY = coordinateLocal[1];

  VLOG(3) << "getNodeNo(" << coordinateLocal << "), nNodesLocalWithoutGhosts: " << this->meshPartition()->nNodesLocalWithoutGhosts()
    << ", nNodesLocalWithoutGhosts: (" << this->meshPartition()->nNodesLocalWithoutGhosts(0) << "," << this->meshPartition()->nNodesLocalWithoutGhosts(1) << ")";

  if (localX == this->meshPartition()->nNodesLocalWithoutGhosts(0))   // point is on right ghost row
  {
    if (localY == this->meshPartition()->nNodesLocalWithoutGhosts(1))  // point is on top ghost row
    {
      // top right ghost point
      VLOG(3) << "  a: " << this->meshPartition()->nNodesLocalWithGhosts()-1;
      return this->meshPartition()->nNodesLocalWithGhosts()-1;
    }
    else
    {
      // on right ghost row
      VLOG(3) << "  b: " << this->meshPartition()->nNodesLocalWithoutGhosts() + localY;
      return this->meshPartition()->nNodesLocalWithoutGhosts() + localY;
    }
  }
  else
  {
    if (localY == this->meshPartition()->nNodesLocalWithoutGhosts(1))   // point is on top ghost row
    {
      // on top ghost row
      if (this->meshPartition()->hasFullNumberOfNodes(0))
      {
        // there are only ghost on the top (y+)
        VLOG(3) << "  c: " << this->meshPartition()->nNodesLocalWithoutGhosts() + localX;
        return this->meshPartition()->nNodesLocalWithoutGhosts() + localX;
      }
      else
      {
        // there are ghosts on the right (x+) and top (y+)
        VLOG(3) << "  d: " << this->meshPartition()->nNodesLocalWithoutGhosts() + this->meshPartition()->nNodesLocalWithoutGhosts(1) + localX;
        return this->meshPartition()->nNodesLocalWithoutGhosts() + this->meshPartition()->nNodesLocalWithoutGhosts(1) + localX;
      }
    }
    else
    {
      // point is in interior
      VLOG(3) << " e: " << this->meshPartition()->nNodesLocalWithoutGhosts() + this->meshPartition()->nNodesLocalWithoutGhosts(0)*localY + localX;
      return this->meshPartition()->nNodesLocalWithoutGhosts(0)*localY + localX;
    }
  }
}

// get local node no from the coordinates (x,y,z), may be a ghost node, for 3D
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>>::
getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const
{
  dof_no_t localX = coordinateLocal[0];
  dof_no_t localY = coordinateLocal[1];
  dof_no_t localZ = coordinateLocal[2];

  VLOG(3) << "getNodeNo(" << coordinateLocal << "), nNodesLocalWithoutGhosts: " << this->meshPartition()->nNodesLocalWithoutGhosts()
    << ", nNodesLocalWithoutGhosts: (" << this->meshPartition()->nNodesLocalWithoutGhosts(0) << "," << this->meshPartition()->nNodesLocalWithoutGhosts(1) << "," << this->meshPartition()->nNodesLocalWithoutGhosts(2) << ")"
    << ", hasFullNumberOfNodes: (" << this->meshPartition()->hasFullNumberOfNodes(0) << "," << this->meshPartition()->hasFullNumberOfNodes(1) << "," << this->meshPartition()->hasFullNumberOfNodes(2) << ")";

  if (localX == this->meshPartition()->nNodesLocalWithoutGhosts(0))   // point is on right ghost row (x+)
  {
    if (this->meshPartition()->hasFullNumberOfNodes(1)) // there is no ghosts at y+
    {
      // there are only ghosts at x+
      return this->meshPartition()->nNodesLocalWithoutGhosts()
          + this->meshPartition()->nNodesLocalWithoutGhosts(1)*localZ + localY;
    }
    else  // there are ghosts at y+
    {
      const int nNodesPerL = this->meshPartition()->nNodesLocalWithoutGhosts(0) + this->meshPartition()->nNodesLocalWithGhosts(1);
      if (localY == this->meshPartition()->nNodesLocalWithoutGhosts(1))  // point is on back ghost row (y+)
      {
        // ghost is at right back vertical edge
        return this->meshPartition()->nNodesLocalWithoutGhosts()
          + nNodesPerL*localZ + this->meshPartition()->nNodesLocalWithoutGhosts(1) + this->meshPartition()->nNodesLocalWithoutGhosts(0);
      }
      else
      {
        // there are ghosts at x+ and y+
        return this->meshPartition()->nNodesLocalWithoutGhosts()
            + nNodesPerL*localZ + localY;
      }
    }
  }
  else if (localY == this->meshPartition()->nNodesLocalWithoutGhosts(1))  // point is on back ghost row (y+)
  {
    if (this->meshPartition()->hasFullNumberOfNodes(0)) // there is no ghosts at x+
    {
      // there are only ghosts at y+
      return this->meshPartition()->nNodesLocalWithoutGhosts()
          + this->meshPartition()->nNodesLocalWithoutGhosts(0)*localZ + localX;
    }
    else  // there are ghosts at x+
    {
      // there are ghosts at x+ and y+
      const int nNodesPerL = this->meshPartition()->nNodesLocalWithoutGhosts(0) + this->meshPartition()->nNodesLocalWithGhosts(1);
      return this->meshPartition()->nNodesLocalWithoutGhosts()
          + nNodesPerL*localZ + this->meshPartition()->nNodesLocalWithoutGhosts(1) + localX;
    }
  }
  else if (localZ == this->meshPartition()->nNodesLocalWithoutGhosts(2))  // point is on top ghost row
  {
    if (this->meshPartition()->hasFullNumberOfNodes(1)) // there is no ghosts at y+
    {
      if (this->meshPartition()->hasFullNumberOfNodes(0)) // there is no ghosts at x+
      {
        // there are only ghosts at z+ and point is at z+
        return this->meshPartition()->nNodesLocalWithoutGhosts()
          + this->meshPartition()->nNodesLocalWithoutGhosts(0)*localY + localX;
      }
      else  // there are ghosts at x+
      {
        // there are ghosts at x+ and z+ and point is at z+
        return this->meshPartition()->nNodesLocalWithoutGhosts()
          + this->meshPartition()->nNodesLocalWithoutGhosts(1)*this->meshPartition()->nNodesLocalWithoutGhosts(2)
          + this->meshPartition()->nNodesLocalWithGhosts(0)*localY + localX;
      }
    }
    else  // there are ghosts at y+
    {
      if (this->meshPartition()->hasFullNumberOfNodes(0)) // there is no ghosts at x+
      {
        // there are ghosts at y+ and z+ and point is at z+
        return this->meshPartition()->nNodesLocalWithoutGhosts()
          + this->meshPartition()->nNodesLocalWithoutGhosts(0)*this->meshPartition()->nNodesLocalWithoutGhosts(2)
          + this->meshPartition()->nNodesLocalWithoutGhosts(0)*localY + localX;
      }
      else  // there are ghosts at x+
      {
        // there are ghosts at x+, y+ and z+ and point is at z+
        const int nNodesPerL = this->meshPartition()->nNodesLocalWithoutGhosts(0) + this->meshPartition()->nNodesLocalWithGhosts(1);
        return this->meshPartition()->nNodesLocalWithoutGhosts()
          + nNodesPerL*this->meshPartition()->nNodesLocalWithoutGhosts(2)
          + this->meshPartition()->nNodesLocalWithGhosts(0)*localY + localX;
      }
    }
  }
  else
  {
    // ghost is in interior
    return this->meshPartition()->nNodesLocalWithoutGhosts(0)*this->meshPartition()->nNodesLocalWithoutGhosts(1)*localZ
      + this->meshPartition()->nNodesLocalWithoutGhosts(0)*localY + localX;
  }
}


};  // namespace
