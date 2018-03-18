#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/02_field_variable_set_get_structured.h"
#include "field_variable/field_variable_data.h"
#include "basis_on_mesh/05_basis_on_mesh.h"

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> :
  public FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! inherited constructor 
  using FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::FieldVariableSetGetStructured;
  
  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo);

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2);

  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
  Vec &values();
  
  //! get the number of components
  int nComponents() const;
  
  //! get the names of the components
  std::vector<std::string> componentNames() const;
  
  //! get the number of elements in the coordinate directions
  std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElementsPerCoordinateDirection() const;
  
  //! get the total number of elements
  element_no_t nElements() const;
};

};  // namespace

#include "field_variable/structured/03_field_variable_data_structured_deformable.tpp"
