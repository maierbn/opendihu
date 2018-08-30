#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/05_field_variable_data_structured_deformable.h"
#include "field_variable/field_variable_set_get.h"
#include "function_space/function_space.h"

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableData<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! inherited constructor
  using FieldVariableData<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableData;

  //! avoid name hiding of "setValues" method
  using FieldVariableData<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::setValues;

  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &rhs);
};

};  // namespace

#include "field_variable/structured/06_field_variable_set_get_structured_deformable.tpp"
