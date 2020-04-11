#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"
#include "field_variable/field_variable_set_get.h"
#include "function_space/function_space.h"

namespace FieldVariable
{

/** FieldVariable class for RegularFixed mesh
 */
template<typename FunctionSpaceType, int nComponents>
class FieldVariableSetGetRegularFixed :
  public FieldVariableSetGetComponent<FunctionSpaceType,nComponents>
{
public:

  //! inherited constructor
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::FieldVariableSetGetComponent;

  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::getValues;
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::getElementValues;
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::getValue;
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::getValuesWithGhosts;
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::getValuesWithoutGhosts;
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::setValue;
  using FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::setValues;

  //! for a specific component, get all values
  void getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false) const;

  //! for a specific component, get all values
  void getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(std::array<std::vector<double>,nComponents> &values, bool onlyNodalValues=false) const;

  //! for a specific component, get values from their local dof no.s
  template<int N>
  void getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const;

  //! for a specific component, get values from their local dof no.s, as vector
  void getValues(int componentNo, const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const;

  //! get values from their local dof no.s for all components
  template<int N>
  void getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values) const;

  //! get values from their local dof no.s for all components
  void getValues(std::vector<dof_no_t> dofLocalNo, std::vector<std::array<double,nComponents>> &values) const;

  //! for a specific component, get the values corresponding to all element-local dofs
  void getElementValues(int componentNo, element_no_t elementNoLocal, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNoLocal, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! get the values corresponding to all element-local dofs for all components, vectorized version for Vc::double_v::size() elements at once
  void getElementValues(Vc::int_v elementNoLocal, std::array<std::array<Vc::double_v,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! for a specific component, get a single value from local dof no.
  double getValue(int componentNo, node_no_t dofLocalNo) const;

  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<FunctionSpaceType,nComponents> &rhs);
private:

  //! get the derivative dof value for the geometry field, component componentNo (x,y or z) when using Hermite polynomials
  double getGeometryFieldHermiteDerivative(int nodeLocalDofIndex, int componentNo) const;
};

} // namespace

#include "field_variable/structured/06_field_variable_set_get_structured_regular_fixed.tpp"
