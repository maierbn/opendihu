#include "field_variable/structured/03_field_variable_data_structured_deformable.h"

#include <sstream>
#include <cassert>

#include "utility/string_utility.h"

namespace FieldVariable
{
 
//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo)
{
  // use the implementation of FieldVariableDataStructured
  FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
    outputHeaderExelem(file, currentElementGlobalNo);
}

//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  // use the implementation of FieldVariableDataStructured
  FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
    outputHeaderExnode(file, currentNodeGlobalNo, valueIndex);
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<int D, typename BasisFunctionType>
bool FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
    haveSameExfileRepresentation(element1, element2);
}

//! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
template<int D, typename BasisFunctionType>
Vec &FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
values()
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::values();
}

//! get the number of components
template<int D, typename BasisFunctionType>
int FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
nComponents() const
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::nComponents();
}

//! get the number of elements
template<int D, typename BasisFunctionType>
std::array<element_no_t, BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::Mesh::dim()> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
nElementsPerCoordinateDirection() const
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::nElementsPerCoordinateDirection();
}

//! get the number of elements
template<int D, typename BasisFunctionType>
element_no_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
nElements() const
{
  // use the implementation of FieldVariableDataStructured
  return this->nElements();
}

//! get the names of the components that are part of this field variable
template<int D, typename BasisFunctionType>
std::vector<std::string> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
componentNames() const
{
  // use the implementation of FieldVariableDataStructured
  return FieldVariableDataStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::componentNames();
}
};