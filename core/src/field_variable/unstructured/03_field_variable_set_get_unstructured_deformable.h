#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/unstructured/02_field_variable_data_unstructured_deformable.h"
#include "field_variable/unstructured/component.h"
#include "field_variable/unstructured/element_to_node_mapping.h"
#include "field_variable/unstructured/node_to_dof_mapping.h"
#include "function_space/function_space.h"
#include "function_space/06_function_space_dofs_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/field_variable_set_get.h"

namespace FieldVariable
{

/** FieldVariable class for UnstructuredDeformable mesh
 */
template<typename FunctionSpaceType, int nComponents>
class FieldVariableSetGetUnstructured :
  public FieldVariableData<FunctionSpaceType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableData<FunctionSpaceType,nComponents>::FieldVariableData;

  //! for a specific component, get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false) const;

  //! for a specific component, get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithGhosts(std::array<std::vector<double>,nComponents> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(std::array<std::vector<double>,nComponents> &values, bool onlyNodalValues=false) const;

  //! for a specific component, get values from their local dof no.s, as array, therefore templated by the number of elements, N, to retrieve
  template<int N>
  void getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const;

  //! for a specific component, get values from their local dof no.s, as vector
  void getValues(int componentNo, const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const;

  //! get values for all components, from their local dof no.s, as contiguous vector in order [comp0, comp0, comp0, ..., comp1, comp1, ...]
  void getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const;

  //! get values from their local dof no.s for all components, this eventually does not get all values if there are multiple versions
  template<int N>
  void getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values) const;

  //! for a specific component, get the values corresponding to all element-local dofs
  void getElementValues(int componentNo, element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! for a specific component, get a single value from local dof no.
  double getValue(int componentNo, node_no_t dofLocalNo) const;

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo) const;

  //! copy the values of a given component to a new single-component field variable
  void extractComponentCopy(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable);

  //! extract the specified component from the field variable by using the raw data array in the given field variable. Afterwards this field variable is invalid and can only be used again after restoreExtractedComponent has been called
  void extractComponentShared(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable);

  //! set the values for the given component from a petsc Vec
  void setValues(int componentNo, Vec petscVector);

  //! set the values for the given component from the other field variable
  void setValues(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> fieldVariable);

  //! set values for a given components for given dofs
  void setValues(int componentNo, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for a given components for given dofs
  template<int N>
  void setValues(int componentNo, const std::array<dof_no_t,N> &dofNosLocal, const std::array<double,N> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<FunctionSpaceType,nComponents> &rhs);

  //! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for dofs with a single component, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for all components for N dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  template<int N>
  void setValues(const std::array<dof_no_t,N> &dofNosLocal, const std::array<std::array<double,nComponents>,N> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof for a given component, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(int componentNo, dof_no_t dofLocalNo, double value, InsertMode petscInsertMode);

  //! set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithoutGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set value for all dofs
  void setValues(double value);

  //! set values for the all component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the all component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithoutGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set value to zero for all dofs
  void zeroEntries();
};

} // namespace

#include "field_variable/unstructured/03_field_variable_set_get_unstructured_deformable.tpp"
