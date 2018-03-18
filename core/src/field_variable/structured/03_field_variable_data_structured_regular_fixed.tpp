#include "field_variable/structured/03_field_variable_data_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

template<int D, typename BasisFunctionType>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
setMeshWidth(double meshWidth)
{
  meshWidth_ = meshWidth;
}

template<int D, typename BasisFunctionType>
double FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
meshWidth() const
{
  return this->meshWidth_;
}

//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  // use the implementation of FieldVariableDataStructured
  FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
    outputHeaderExelem(file, currentElementGlobalNo, fieldVariableNo);
}

//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
  // use the implementation of FieldVariableDataStructured
  FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
    outputHeaderExnode(file, currentNodeGlobalNo, valueIndex, fieldVariableNo);
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<int D, typename BasisFunctionType>
bool FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
    haveSameExfileRepresentation(element1, element2);
}

//! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
template<int D, typename BasisFunctionType>
Vec &FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
values()
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
    values();
}

//! get the number of components
template<int D, typename BasisFunctionType>
int FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
nComponents() const
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
    nComponents();
}
/*
//! get the number of elements
template<int D, typename BasisFunctionType>
std::array<element_no_t, BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::Mesh::dim()> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
nElementsPerCoordinateDirection() const
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::nElementsPerCoordinateDirection();
}*/

//! get the number of elements
template<int D, typename BasisFunctionType>
element_no_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
nElements() const
{
  return this->nElements();
}

//! get the names of the components that are part of this field variable
template<int D, typename BasisFunctionType>
std::vector<std::string> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
componentNames() const
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
    componentNames();
}
};
