#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"
#include "field_variable/field_variable_data.h"
#include "function_space/function_space.h"

namespace FieldVariable
{

/** FieldVariable class for RegularFixed mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableSetGetComponent<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! inherited constructor
  using FieldVariableSetGetComponent<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableSetGetComponent;

protected:
};

};  // namespace

#include "field_variable/structured/05_field_variable_data_structured_regular_fixed.tpp"
