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
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableSetGetUnstructured<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:

  //! inherited constructors
  using FieldVariableSetGetUnstructured<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableSetGetUnstructured;

  //! avoid name hiding of "getValue" in parent classes
  using FieldVariableSetGetUnstructured<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::getValue;

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo);


};

/** For 1 component
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1> :
  public FieldVariableSetGetUnstructured<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>
{
public:

  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;

  //! inherited constructors
  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::FieldVariableSetGetUnstructured;

  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::getElementValues;
  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::getValue;
  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::setValue;
  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::setValues;
  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::setValuesWithGhosts;
  using FieldVariableSetGetUnstructured<BasisOnMeshType,1>::setValuesWithoutGhosts;

  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);

  //! get a single value from local dof no. for all components
  double getValue(node_no_t dofLocalNo);

  //! get all stored local values
  void getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues=false);
  
  //! get all stored local values
  void getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues=false);
  
  //! set a single dof (all components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for given dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

};

};  // namespace

#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.tpp"
