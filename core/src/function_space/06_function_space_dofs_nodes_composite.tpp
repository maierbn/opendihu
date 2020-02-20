#include "function_space/06_function_space_dofs_nodes.h"

#include <Python.h>  // has to be the first included header
#include "easylogging++.h"

#include <cmath>
#include <array>
#include <string>

namespace FunctionSpace
{

template<int D,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nNodesLocalWithGhosts() const
{
  return this->meshPartition()->nNodesLocalWithGhosts();
}

template<int D,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nNodesLocalWithoutGhosts() const
{
  return this->meshPartition()->nNodesLocalWithoutGhosts();
}

template<int D,typename BasisFunctionType>
dof_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nDofsLocalWithGhosts() const
{
  return this->meshPartition()->nDofsLocalWithGhosts();
}

template<int D,typename BasisFunctionType>
dof_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nDofsLocalWithoutGhosts() const
{
  return this->meshPartition()->nDofsLocalWithoutGhosts();
}

template<int D,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nNodesGlobal() const
{
  return this->meshPartition()->nNodesGlobal();
}

template<int D,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nDofsGlobal() const
{
  return this->meshPartition()->nDofsGlobal();
}

template<int D,typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  nodes.resize(this->nNodesLocalWithGhosts()*3);

  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nNodesLocalWithGhosts(); nodeGlobalNo++)
  {
    node_no_t nodeFirstDofGlobalNo = this->geometryField_->nodeToDofMapping()->getNodeDofs(nodeGlobalNo)[0];
    Vec3 position = this->geometryField_->template getValue(nodeFirstDofGlobalNo);
    int index = nodeGlobalNo*3;
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}


} // namespace
