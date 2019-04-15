#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "field_variable/structured/03_field_variable_set_get_structured.h"

namespace FieldVariable
{

/** Class that implements get and set methods that are different for different values of nComponent.
 *  For 1 component they use double instead of std::array<double,1>.
 *  For multiple components they use std::array<double,nComponent>.
 */
/* For >1 components
 */
template<typename FunctionSpaceType, int nComponents>
class FieldVariableSetGetComponent :
  public FieldVariableSetGetStructured<FunctionSpaceType,nComponents>
{
public:

  //! inherited constructors
  using FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::FieldVariableSetGetStructured;

  //! avoid name hiding of "getValue" in parent classes
  using FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::getValue;

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo) const;


};

/** For 1 component
 */
template<typename FunctionSpaceType>
class FieldVariableSetGetComponent<FunctionSpaceType,1> :
  public FieldVariableSetGetStructured<FunctionSpaceType,1>
{
public:

  //! inherited constructors
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::FieldVariableSetGetStructured;

  using FieldVariableSetGetStructured<FunctionSpaceType,1>::getElementValues;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::getValue;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::getValues;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::getValuesWithGhosts;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::getValuesWithoutGhosts;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::setValuesWithGhosts;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::setValuesWithoutGhosts;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::setValue;
  using FieldVariableSetGetStructured<FunctionSpaceType,1>::setValues;

  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNoLocal, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! get a single value from local dof no. for all components
  double getValue(node_no_t dofLocalNo) const;

  //! get values from their local dof no.s, as vector
  void getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const;

  //! get all stored local values
  void getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues=false) const;
  
  //! get all stored local values
  void getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues=false) const;
  
  //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set the values from a petsc Vec
  void setValues(Vec petscVector);

  //! set values for the single component for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  template<int nValues>
  void setValues(const std::array<dof_no_t,nValues> dofNosLocal, std::array<double,nValues> values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);
  
  //! set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);
};

}  // namespace

#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.tpp"
