#include "field_variable/unstructured/03_field_variable_set_get_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>
#include <iomanip>
#include <cassert>

namespace FieldVariable
{

using namespace StringUtility;

//! get values from their local dof no.s for all components, this eventually does not get all values if there are multiple versions
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values) const
{
  assert(this->values_);
  
  std::array<double,N*nComponents> resultVector;
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // get values and assign them to result values vector
    this->values_->getValues(componentIndex, N, dofLocalNo.data(), resultVector.data() + componentIndex*N);
  }
  
  // transform local dof no.s to vector indices of first component
  for (int dofIndex = 0; dofIndex < N; dofIndex++)
  {
    // copy retrieved values to output array
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      values[dofIndex][componentIndex] = resultVector[componentIndex*N + dofIndex];
    }
  }
}

//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getValues(values, onlyNodalValues);
}

//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues) const
{
  this->getValuesWithGhosts(componentNo, values, onlyNodalValues);
}

//! for a specific component, get values from their local dof no.s
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].template getValues<N>(dofLocalNo, values);
}

//! for a specific component, get values from their local dof no.s
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofLocalNo, std::vector<double> &values) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getValues(dofLocalNo, values);
}

//! for a specific component, get a single value from local dof no.
template<typename FunctionSpaceType, int nComponents>
double FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getValue(int componentNo, node_no_t dofLocalNo) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  return this->component_[componentNo].getValue(dofLocalNo);
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getElementValues(elementNo, values);
}

//! get the values corresponding to all element-local dofs for all components
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(this->values_);
  
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  const std::vector<dof_no_t> &dofLocalNo = this->elementToDofMapping_->getElementDofs(elementNo);
  
  std::array<double,nDofsPerElement*nComponents> resultVector;
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // get values and assign them to result values vector
    this->values_->getValues(componentIndex, nDofsPerElement, dofLocalNo.data(), resultVector.data() + componentIndex*nDofsPerElement);
  }
  
  // transform local dof no.s to vector indices of first component
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    // copy retrieved values to output array
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      values[dofIndex][componentIndex] = resultVector[componentIndex*nDofsPerElement + dofIndex];
    }
  }
}

//! copy the values from another field variable of the same type
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValues(FieldVariable<FunctionSpaceType,nComponents> &rhs)
{
  this->values_->setValues(*rhs.partitionedPetscVec());
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->nComponents == 1 && nComponents == 1);
  const int nValues = values.size();
  assert(this->values_);

  this->values_->setValues(0, nValues, (const int *) dofNosLocal.data(), values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set value for all dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValues(double value)
{
  // get number of dofs
  assert(this->functionSpace_);
  const dof_no_t nDofs = this->functionSpace_->nDofsLocalWithoutGhosts();

  std::vector<double> valueBuffer(nDofs,value);
  
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->setValuesWithoutGhosts(componentIndex, valueBuffer, INSERT_VALUES);
  }
}

//! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(dofNosLocal.size() == values.size());
  assert(this->values_);
 
  const int nValues = values.size();
  std::vector<double> valuesBuffer(nValues);

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // loop over dofs and prepare values of current component
    for (int dofIndex = 0; dofIndex < nValues; dofIndex++)
    {
      valuesBuffer[dofIndex] = values[dofIndex][componentIndex];
    }
    
    this->values_->setValues(componentIndex, nValues, dofNosLocal.data(), valuesBuffer.data(), petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode)
{
  assert(this->values_);
  
  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->setValues(componentIndex, 1, &dofLocalNo, value.data()+componentIndex, petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValuesWithGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());
  assert(this->values_);
 
  this->values_->setValues(componentNo, values.size(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValuesWithoutGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
  assert(this->values_);
   
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->values_->setValues(componentNo, values.size(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

//! set values for all components for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValuesWithGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());
  
  this->setValues(this->functionSpace_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

//! set values for all components for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
setValuesWithoutGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
  
  this->setValues(this->functionSpace_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

//! set value to zero for all dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,nComponents>::
zeroEntries()
{
  assert(this->values_);
  this->values_->zeroEntries();
}

};  // namespace
