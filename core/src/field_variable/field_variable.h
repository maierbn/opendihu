#pragma once

#include "field_variable/structured/06_field_variable_set_get_structured_deformable.h"
#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.h"
#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.h"

namespace FieldVariable
{

/** General field variable with != 1 components. A field variable is defined on a BasisOnMesh, i.e. knows mesh type and basis function type.
 */
template<typename BasisOnMeshType,int nComponents>
class FieldVariable :
  public FieldVariableSetGet<BasisOnMeshType,nComponents>
{ 
public:
  //! inherited constructor 
  using FieldVariableSetGet<BasisOnMeshType,nComponents>::FieldVariableSetGet;
  
  typedef BasisOnMeshType BasisOnMesh;
};

/** General scalar field variable.
 * A field variable is defined on a BasisOnMesh, i.e. knows mesh type and basis function type.
 * Scalar field variables can compute a gradient field.
 */
template<typename BasisOnMeshType>
class FieldVariable<BasisOnMeshType,1> : 
  public FieldVariableSetGet<BasisOnMeshType,1>
{
public:
  //! inherited constructor 
  using FieldVariableSetGet<BasisOnMeshType,1>::FieldVariableSetGet;
  
  typedef BasisOnMeshType BasisOnMesh;
 
  //! fill the gradient field with the gradient values in world coordinates of this field variable. This is only possible for scalar fields.
  void computeGradientField(FieldVariable<BasisOnMeshType, BasisOnMeshType::dim()> &gradientField);
  
};

// output operator
template<typename BasisOnMeshType,int nComponents>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  rhs.output(stream);
  return stream;
}

};  // namespace

#include "field_variable/07_field_variable_gradient.tpp"