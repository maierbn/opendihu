#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/06_field_variable_set_get_structured_regular_fixed.h"
#include "field_variable/field_variable_set_get.h"
#include "function_space/function_space.h"

namespace FieldVariable
{

/** FieldVariable class for RegularFixed mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableSetGetRegularFixed<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! inherited constructor
  using FieldVariableSetGetRegularFixed<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableSetGetRegularFixed;
  using FieldVariableSetGetRegularFixed<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::getValue;

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo) const;
};

/** Partial specialization for single component field variable
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,1> :
  public FieldVariableSetGetRegularFixed<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,1>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! inherited constructor
  using FieldVariableSetGetRegularFixed<FunctionSpaceType,1>::FieldVariableSetGetRegularFixed;
  
};

} // namespace

#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.tpp"
