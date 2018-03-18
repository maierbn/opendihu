#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/03_field_variable_data_structured_deformable.h"
#include "field_variable/field_variable_set_get.h"
#include "basis_on_mesh/05_basis_on_mesh.h"

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> :
  public FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! inherited constructor 
  using FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::FieldVariableData;
  
  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &rhs);
};

};  // namespace

#include "field_variable/structured/04_field_variable_set_get_structured_deformable.tpp"
