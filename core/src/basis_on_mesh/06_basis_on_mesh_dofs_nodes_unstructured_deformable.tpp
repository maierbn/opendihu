#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes.h"

#include <Python.h>  // has to be the first included header
#include "easylogging++.h"

#include <cmath>
#include <array>
#include <string>

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesLocalWithGhosts() const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return this->geometryField_->nNodes();
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesLocalWithoutGhosts() const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return this->geometryField_->nNodes();
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofsLocalWithGhosts() const
{
  // parallelism not implemented for unstructured grids, therefore there is no distinction between global and local numbers
  return this->nDofs_;
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofsLocalWithoutGhosts() const
{
  // parallelism not implemented for unstructured grids, therefore there is no distinction between global and local numbers
  return this->nDofs_;
}

template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodesGlobal() const
{
  // parallelism not implemented for unstructured grids, therefore there is no distinction between global and local numbers
  assert(this->geometryField_);
  return this->geometryField_->nNodes();
}

template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofsGlobal() const
{
  // parallelism not implemented for unstructured grids
  return this->nDofs_;
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
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


};  // namespace
