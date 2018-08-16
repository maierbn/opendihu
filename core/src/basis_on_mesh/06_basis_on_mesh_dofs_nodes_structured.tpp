#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes_structured.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

namespace BasisOnMesh
{
  
template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodesStructuredStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesLocalWithGhosts(int coordinateDirection) const
{
  assert(this->meshPartition_->localSize(coordinateDirection) == this->nElementsPerCoordinateDirectionLocal(coordinateDirection));
 
  return this->meshPartition_->nNodesLocalWithGhosts(coordinateDirection);
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodesStructuredStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesLocalWithoutGhosts(int coordinateDirection) const
{
  assert(this->meshPartition_->localSize(coordinateDirection) == this->nElementsPerCoordinateDirectionLocal(coordinateDirection));
 
  return this->meshPartition_->nNodesLocalWithoutGhosts(coordinateDirection);
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodesStructuredStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesLocalWithGhosts() const
{
  return this->meshPartition_->nNodesLocalWithGhosts();
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodesStructuredStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesLocalWithoutGhosts() const
{
  return this->meshPartition_->nNodesLocalWithoutGhosts();
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofsNodesStructuredStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofsLocalWithGhosts() const
{
  return nNodesLocalWithGhosts() * this->nDofsPerNode();
}

template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesGlobal() const
{
  this->meshPartition_->nNodesGlobal();
}

template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesGlobal(int coordinateDirection) const
{
  return this->meshPartition_->nNodesGlobal(coordinateDirection);
}

template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofsGlobal() const
{
  return nNodesGlobal() * this->nDofsPerNode();
}

//! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  nodes.resize(this->nDofsLocalWithGhosts()*3);

  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nDofsLocalWithGhosts(); nodeGlobalNo++)
  {

    node_no_t firstNodeDofGlobalNo = nodeGlobalNo*this->nDofsPerNode();

    int index = nodeGlobalNo*3;
    Vec3 position = this->geometryField_->getValue(firstNodeDofGlobalNo);
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}


};  // namespace