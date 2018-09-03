#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/unstructured/03_field_variable_set_get_unstructured_deformable.h"

namespace FieldVariable
{

/** Class that implements get and set methods that are different for different values of nComponent.
 *  For 1 component they use double instead of std::array<double,1>.
 *  For multiple components they use std::array<double,nComponent>.
 */
/* For >1 components
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableSetGetUnstructured<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:

  //! inherited constructors
  using FieldVariableSetGetUnstructured<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableSetGetUnstructured;

  //! avoid name hiding of "getValue" in parent classes
  using FieldVariableSetGetUnstructured<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::getValue;

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo) const;

};

/** For 1 component
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1> :
  public FieldVariableSetGetUnstructured<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>
{
public:

  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! inherited constructors
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::FieldVariableSetGetUnstructured;

  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::getElementValues;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::getValue;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::getValuesWithGhosts;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::getValuesWithoutGhosts;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::setValue;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::setValues;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::setValuesWithGhosts;
  using FieldVariableSetGetUnstructured<FunctionSpaceType,1>::setValuesWithoutGhosts;

  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const;

  //! get a single value from local dof no. for all components
  double getValue(node_no_t dofLocalNo) const;

  //! get all stored local values
  void getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues=false) const;
  
  //! get all stored local values
  void getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues=false) const;
  
  //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for given dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

};

};  // namespace

#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.tpp"
