#pragma once

#include "field_variable/structured/06_field_variable_set_get_structured_deformable.h"
#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.h"
#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.h"
#include "mesh/type_traits.h"

namespace FieldVariable
{

/** General field variable with != 1 components. A field variable is defined on a FunctionSpace, i.e. knows mesh type and basis function type.
 */
template<typename FunctionSpaceType,int nComponents,typename Mesh=typename FunctionSpaceType::Mesh>
class FieldVariableGradient :
  public FieldVariableSetGet<FunctionSpaceType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableSetGet<FunctionSpaceType,nComponents>::FieldVariableSetGet;
};

/** Scalar field variable with structured mesh.
 * A field variable is defined on a FunctionSpace, i.e. knows mesh type and basis function type.
 * Scalar field variables can compute a gradient field.
 */
template<typename FunctionSpaceType>
class FieldVariableGradient<FunctionSpaceType,1,::Mesh::isStructured<typename FunctionSpaceType::Mesh>> :
  public FieldVariableSetGet<FunctionSpaceType,1>
{
public:
  //! inherited constructor
  using FieldVariableSetGet<FunctionSpaceType,1>::FieldVariableSetGet;

  //! fill the gradient field with the gradient values in world coordinates of this field variable. This is only possible for scalar fields.
  void computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField,
                            std::shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber = nullptr);

};

/** Scalar field variable with unstructured mesh.
 * A field variable is defined on a FunctionSpace, i.e. knows mesh type and basis function type.
 * Scalar field variables can compute a gradient field.
 */
template<typename FunctionSpaceType>
class FieldVariableGradient<FunctionSpaceType,1,::Mesh::UnstructuredDeformableOfDimension<FunctionSpaceType::dim()>> :
  public FieldVariableSetGet<FunctionSpaceType,1>
{
public:
  //! inherited constructor
  using FieldVariableSetGet<FunctionSpaceType,1>::FieldVariableSetGet;

  //! fill the gradient field with the gradient values in world coordinates of this field variable. This is only possible for scalar fields.
  void computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField,
                            std::shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber = nullptr);

};

} // namespace

#include "field_variable/08_field_variable_gradient.tpp"
