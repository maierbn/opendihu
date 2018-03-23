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
template<typename BasisOnMeshType, int nComponents>
class FieldVariableSetGetComponent : 
  public FieldVariableSetGetStructured<BasisOnMeshType,nComponents>
{
public:

  //! inherited constructors
  using FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::FieldVariableSetGetStructured;
 
  //! avoid name hiding of "getValue" in parent classes
  using FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::getValue;
  
  //! get a single value from global dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);


};

/** For 1 component
 */
template<typename BasisOnMeshType>
class FieldVariableSetGetComponent<BasisOnMeshType,1> : 
  public FieldVariableSetGetStructured<BasisOnMeshType,1>
{
public:
  
  //! inherited constructors
  using FieldVariableSetGetStructured<BasisOnMeshType,1>::FieldVariableSetGetStructured;
 
  using FieldVariableSetGetStructured<BasisOnMeshType,1>::getElementValues;
  using FieldVariableSetGetStructured<BasisOnMeshType,1>::getValue;
  using FieldVariableSetGetStructured<BasisOnMeshType,1>::setValue;
  using FieldVariableSetGetStructured<BasisOnMeshType,1>::setValues;
  
  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get a single value from global dof no. for all components
  double getValue(node_no_t dofGlobalNo);
  
  //! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  void setValue(dof_no_t dofGlobalNo, double value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

};

};   // namespace

#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.tpp"
