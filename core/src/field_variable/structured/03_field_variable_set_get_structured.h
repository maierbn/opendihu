#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>
#include <Vc/Vc>

#include "field_variable/structured/02b_field_variable_data_structured_for_surface.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

namespace FieldVariable
{

/** Field variable for a structured mesh, i.e. dof and node information are purely implicit.
 *  This is used for RegularFixed and StructuredDeformable meshes.
 */
template<typename FunctionSpaceType, int nComponents>
class FieldVariableSetGetStructured :
  public FieldVariableDataStructuredForSurface<FunctionSpaceType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableDataStructuredForSurface<FunctionSpaceType,nComponents>::FieldVariableDataStructuredForSurface;

  using FieldVariableDataStructuredForSurface<FunctionSpaceType,nComponents>::setValues;

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
  void getValuesWithoutGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues=false) const;

  //! get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(std::array<std::vector<double>,nComponents> &values, bool onlyNodalValues=false) const;

  //! for a specific component, get values from their local dof no.s, as array, therefore templated by the number of elements, N, to retrieve
  template<int N>
  void getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const;

  //! for a specific component, get values from their local dof no.s, as vector
  void getValues(int componentNo, const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const;

  //! for a specific component, get values from their local dof no.s
  void getValues(int componentNo, int nValues, const dof_no_t *dofLocalNo, std::vector<double> &values) const;

  //! get values for all components, from their local dof no.s, as contiguous vector in order [comp0, comp0, comp0, ..., comp1, comp1, ...]
  void getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const;

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

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo) const;

  //! extract the specified component from the field variable (by copying it) and store it in the given field variable (which already has the data allocated)
  void extractComponentCopy(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable);

  //! extract the specified component from the field variable by using the raw data array in the given field variable. Afterwards this field variable is invalid and can only be used again after restoreExtractedComponent has been called
  void extractComponentShared(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable);

  //! if it is possible to extract a component in the shared way
  bool isExtractComponentSharedPossible(int componentNo);

  //! restore the extracted raw array to petsc and make the field variable usable again, the restored array was from componentNo
  template<int nComponents2>
  void restoreExtractedComponent(std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents2>> extractedVec, int componentNo);

  //! set the values for the given component from a petsc Vec
  void setValues(int componentNo, Vec petscVector);

  //! set the values for the given component from the other field variable
  void setValues(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> fieldVariable);

  //! set values for a given component for given dofs
  void setValues(int componentNo, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for a given component for given dofs
  template<int N>
  void setValues(int componentNo, const std::array<dof_no_t,N> &dofNosLocal, const std::array<double,N> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for a given component for given dofs, using raw pointers
  void setValues(int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for all components for N dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  template<int N>
  void setValues(const std::array<dof_no_t,N> &dofNosLocal, const std::array<std::array<double,nComponents>,N> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for all components for dofs, only nValues values will be set despite potentially more dofNosLocal, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(int nValues, const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof (all components), after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof (all components), after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(Vc::int_v dofLocalNo, const std::array<Vc::double_v,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof for a given component, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(int componentNo, dof_no_t dofLocalNo, double value, InsertMode petscInsertMode);

  //! set a given component of Vc::double_v::size() dofs with the vectorized value, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(int componentNo, Vc::int_v dofLocalNo, Vc::double_v value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a given component of Vc::double_v::size() dofs with the same value
  void setValue(int componentNo, Vc::int_v dofLocalNo, double value, InsertMode petscInsertMode=INSERT_VALUES);

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

  //! set values for the all component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValuesWithoutGhosts(const std::array<std::vector<double>,nComponents> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set value to zero for all dofs
  void zeroEntries();
};

} // namespace

#include "field_variable/structured/03_field_variable_set_get_structured.tpp"
