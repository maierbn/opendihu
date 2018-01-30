#include "basis_on_mesh/04_basis_on_mesh_nodes.h"

#include "easylogging++.h"

#include <cmath>
#include <array>
#include <string>

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
node_idx_t BasisOnMeshNodes<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
nNodes() const
{
  // assert that geometry field variable is set
  assert (this->fieldVariable_.find("geometry") != this->fieldVariable_.end());
  
  return this->fieldVariable_.at("geometry")->nNodes();
}

template<int D,typename BasisFunctionType>
Vec3 BasisOnMeshNodes<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getGeometry(node_idx_t dofGlobalNo) const
{
  // assert that geometry field variable is set
  assert (this->fieldVariable_.find("geometry") != this->fieldVariable_.end());
  
  Vec3 result = this->fieldVariable_.at("geometry")->template getValue<3>(dofGlobalNo);
  return result;
}  
  
//! return an array containing all geometry entries for an element
template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getElementGeometry(element_idx_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement()> &values)
{
  // assert that geometry field variable is set
  assert (this->fieldVariable_.find("geometry") != this->fieldVariable_.end());
  
  const int nDofsPerElement = BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement();
  this->fieldVariable_.at("geometry")->template getElementValues<nDofsPerElement,3>(elementNo, values);
}

//! create a non-geometry field field variable with no values being set, with given component names
template<int D,typename BasisFunctionType>
typename BasisOnMeshNodes<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::FieldVariableType &BasisOnMeshNodes<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
geometryField()
{
  // assert that geometry field variable is set
  assert (this->fieldVariable_.find("geometry") != this->fieldVariable_.end());
  
  return *this->fieldVariable_["geometry"];
}

template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  // assert that geometry field variable is set
  assert (this->fieldVariable_.find("geometry") != this->fieldVariable_.end());
  
  nodes.resize(this->nNodes()*3);
 
  for (node_idx_t nodeGlobalNo = 0; nodeGlobalNo < this->nNodes(); nodeGlobalNo++)
  {
    int nodeFirstDofGlobalNo = this->fieldVariable_.at("geometry")->nodeToDofMapping()->getNodeDofs(nodeGlobalNo)[0];
    Vec3 position = this->fieldVariable_.at("geometry")->template getValue<3>(nodeFirstDofGlobalNo);
    int index = nodeGlobalNo*3;
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}

};  // namespace
