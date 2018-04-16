#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/05_field_variable_data_structured_deformable.h"
#include "field_variable/field_variable_set_get.h"
#include "basis_on_mesh/basis_on_mesh.h"

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! inherited constructor 
  using FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableData;
  
  //! avoid name hiding of "setValues" method
  using FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::setValues;
  
  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &rhs);
};

};  // namespace

#include "field_variable/structured/06_field_variable_set_get_structured_deformable.tpp"
