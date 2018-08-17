#include "basis_on_mesh/05_basis_on_mesh_geometry.h"

#include "easylogging++.h"

namespace BasisOnMesh
{

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
Vec3 BasisOnMeshGeometry<MeshType,BasisFunctionType,DummyForTraits>::
getGeometry(node_no_t dofGlobalNo) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return this->geometryField_->template getValue(dofGlobalNo);
}

//! return an array containing all geometry entries for an element
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
void BasisOnMeshGeometry<MeshType,BasisFunctionType,DummyForTraits>::
getElementGeometry(element_no_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &values)
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  assert (elementNo >= 0);
  assert (elementNo < this->nElementsLocal());

  this->geometryField_->getElementValues(elementNo, values);
}

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
bool BasisOnMeshGeometry<MeshType,BasisFunctionType,DummyForTraits>::
hasGeometryField()
{
  return this->geometryField_ != nullptr;
}

//! create a non-geometry field field variable with no values being set, with given component names
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
typename BasisOnMeshGeometry<MeshType,BasisFunctionType,DummyForTraits>::GeometryFieldType &
BasisOnMeshGeometry<MeshType,BasisFunctionType,DummyForTraits>::
geometryField()
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return *this->geometryField_;
}

};  // namespace
