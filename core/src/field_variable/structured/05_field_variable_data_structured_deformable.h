#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"
#include "field_variable/field_variable_data.h"
#include "basis_on_mesh/basis_on_mesh.h"

namespace BasisOnMesh
{

};

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableData<::BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableSetGetComponent<::BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;

  //! inherited constructor
  using FieldVariableSetGetComponent<::BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableSetGetComponent;

};

};  // namespace

#include "field_variable/structured/05_field_variable_data_structured_deformable.tpp"
