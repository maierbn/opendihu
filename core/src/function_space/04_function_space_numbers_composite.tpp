#include "function_space/04_function_space_numbers_composite.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

#include "mesh/composite.h"

namespace FunctionSpace
{

// element-local dofIndex to local dofNo
template<int D,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getDofNo(element_no_t elementNoLocal, int dofIndex) const
{
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;
  
  return getNodeNo(elementNoLocal, nodeIndex)*nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node
template<int D,typename BasisFunctionType>
void FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + FunctionSpaceBaseDim<D,BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<D,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node
template<int D,typename BasisFunctionType>
void FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<D,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int dofIndex = 0; dofIndex < FunctionSpaceBaseDim<D,BasisFunctionType>::nDofsPerNode(); dofIndex++)
  {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node
template<int D,typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return FunctionSpaceBaseDim<D,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local nodeIndex to local nodeNo
template<int D,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodeNo(element_no_t elementNoLocal, int nodeIndex) const
{
  int subMeshNo = 0;
  element_no_t elementOnMeshNoLocal = 0;
  this->meshPartition_->getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshNo, elementOnMeshNoLocal);

  return this->subFunctionSpaces_[subMeshNo]->getNodeNo(elementOnMeshNoLocal, nodeIndex);
}

// local node no of neighbour node, may be a ghost node,
template<int D,typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const
{
  assert (nodeNoLocal < this->nNodesLocalWithoutGhosts());

  // get all submesh + nodeNo combinations of the current node, since the node could be shared by multiple submeshes
  std::vector<std::pair<int,node_no_t>> subMeshesWithNodes;
  this->meshPartition_->getSubMeshesWithNodes(nodeNoLocal, subMeshesWithNodes);

  // loop over submesh + nodeNos
  for (const std::pair<int,node_no_t>& subMeshWithNodes : subMeshesWithNodes)
  {
    int subMeshNo = subMeshWithNodes.first;
    node_no_t subMeshNodeNo = subMeshWithNodes.second;

    // check if there is a neighbouring node in the current submeseh
    int neighbourNodeNoLocal = this->subFunctionSpaces_[subMeshNo]->getNeighbourNodeNoLocal(subMeshNodeNo, direction);

    // if yes, transform the node no to the composite numbering
    if (neighbourNodeNoLocal != -1)
    {
      return this->meshPartition_->getNodeNoLocalFromSubmesh(subMeshNodeNo, subMeshNodeNo);
    }
  }

  // no neighbour node found
  return -1;
}

template<int D,typename BasisFunctionType>
global_no_t FunctionSpaceNumbersCommon<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodeNoGlobalNatural(element_no_t elementNoLocal, int nodeIndex) const
{
  // this is implemented in the mesh partition
  return this->meshPartition_->getNodeNoGlobalNatural(elementNoLocal, nodeIndex);
}

} // namespace
