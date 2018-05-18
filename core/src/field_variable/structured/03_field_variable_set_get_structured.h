#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "field_variable/structured/02_field_variable_data_structured.h"

namespace FieldVariable
{

/** Field variable for a structured mesh, i.e. dof and node information are purely implicit.
 *  This is used for RegularFixed and StructuredDeformable meshes.
 */
template<typename BasisOnMeshType, int nComponents>
class FieldVariableSetGetStructured :
  public FieldVariableDataStructured<BasisOnMeshType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableDataStructured<BasisOnMeshType,nComponents>::FieldVariableDataStructured;

  //! for a specific component, get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValues(int componentNo, std::vector<double> &values, bool onlyNodalValues=false);

  //! for a specific component, get values from their global dof no.s, as array, therefore templated by the number of elements, N, to retrieve
  template<int N>
  void getValues(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);

  //! for a specific component, get values from their global dof no.s, as vector
  void getValues(int componentNo, std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values);

  //! get values from their global dof no.s for all components
  template<int N>
  void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values);

  //! for a specific component, get the values corresponding to all element-local dofs
  void getElementValues(int componentNo, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);

  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values);

  //! for a specific component, get a single value from global dof no.
  double getValue(int componentNo, node_no_t dofGlobalNo);

  //! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set value for all dofs
  void setValues(double value);

  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  void flushSetValues();
};

};  // namespace

#include "field_variable/structured/03_field_variable_set_get_structured.tpp"
