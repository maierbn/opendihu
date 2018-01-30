#include "field_variable/field_variable_structured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>

namespace FieldVariable
{
 
/*
//! for a specific component, get a single value from global dof no.
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
getValue(std::string component, node_no_t dofGlobalNo)
{
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::getValue(component, dofGlobalNo);
}
*/
 
//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo)
{
  // use the implementation of FieldVariableStructured
  FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
    outputHeaderExelem(file, currentElementGlobalNo);
}

//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  // use the implementation of FieldVariableStructured
  FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
    outputHeaderExnode(file, currentNodeGlobalNo, valueIndex);
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<int D, typename BasisFunctionType>
bool FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
    haveSameExfileRepresentation(element1, element2);
}

//! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
template<int D, typename BasisFunctionType>
Vec &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
values()
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::values();
}

//! get the number of components
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
nComponents() const
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::nComponents();
}

//! get the number of elements
template<int D, typename BasisFunctionType>
std::array<element_no_t, BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>::Mesh::dim()>  FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
nElementsPerDimension() const
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::nElementsPerDimension();
}

//! get the number of elements
template<int D, typename BasisFunctionType>
element_no_t FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
nElements() const
{
  // use the implementation of FieldVariableStructured
  return this->nElements();
}

//! get a single value from global dof no. for all components
template<int D, typename BasisFunctionType>
template<int nComponents>
std::array<double,nComponents> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
getValue(node_no_t dofGlobalNo)
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMeshType>::template getValue<nComponents>(dofGlobalNo);
}


//! get the names of the components that are part of this field variable
template<int D, typename BasisFunctionType>
std::vector<std::string> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
componentNames() const
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::componentNames();
}
};