#include "function_space/05_function_space_geometry.h"

#include "easylogging++.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
Vec3 FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
getGeometry(node_no_t dofGlobalNo) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return this->geometryField_->getValue(dofGlobalNo);
}

//! return an array containing all geometry entries for an element
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
void FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
getElementGeometry(element_no_t elementNo, std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &values)
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  assert (elementNo >= 0);
  if (elementNo >= this->nElementsLocal())
    LOG(ERROR) << "FunctionSpace::getElementGeometry elementNo: " << elementNo << ", nElementsLocal: " << this->nElementsLocal();
  assert (elementNo < this->nElementsLocal());

  this->geometryField_->getElementValues(elementNo, values);
}

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
bool FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
hasGeometryField()
{
  return this->geometryField_ != nullptr;
}

//! create a non-geometry field field variable with no values being set, with given component names
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
typename FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::GeometryFieldType &
FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
geometryField()
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return *this->geometryField_;
}

} // namespace
