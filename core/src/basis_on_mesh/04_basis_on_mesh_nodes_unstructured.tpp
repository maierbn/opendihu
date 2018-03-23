#include "basis_on_mesh/04_basis_on_mesh_nodes.h"

#include <Python.h>  // has to be the first included header
#include "easylogging++.h"

#include <cmath>
#include <array>
#include <string>

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodes() const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  
  return this->geometryField_->nNodes();
}

template<int D,typename BasisFunctionType>
Vec3 BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getGeometry(node_no_t dofGlobalNo) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  
  Vec3 result = this->geometryField_->template getValue(dofGlobalNo);
  return result;
}  
  
//! return an array containing all geometry entries for an element
template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getElementGeometry(element_no_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement()> &values)
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  
  this->geometryField_->getElementValues(elementNo, values);
}


template<int D,typename BasisFunctionType>
bool BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
hasGeometryField()
{
  return this->geometryField_ != nullptr;
}

//! create a non-geometry field field variable with no values being set, with given component names
template<int D,typename BasisFunctionType>
typename BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::GeometryFieldType &BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
geometryField()
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  
  return *this->geometryField_;
}

template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  
  nodes.resize(this->nNodes()*3);
 
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nNodes(); nodeGlobalNo++)
  {
    int nodeFirstDofGlobalNo = this->geometryField_->nodeToDofMapping()->getNodeDofs(nodeGlobalNo)[0];
    Vec3 position = this->geometryField_->template getValue(nodeFirstDofGlobalNo);
    int index = nodeGlobalNo*3;
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}

};  // namespace
