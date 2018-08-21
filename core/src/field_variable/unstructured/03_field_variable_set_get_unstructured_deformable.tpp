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
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values)
{
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
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getValues(values, onlyNodalValues);
}

//! for a specific component, get all values
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  this->getValuesWithGhosts(componentNo, values, onlyNodalValues);
}

//! for a specific component, get values from their local dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].template getValues<N>(dofLocalNo, values);
}

//! for a specific component, get values from their local dof no.s
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofLocalNo, std::vector<double> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getValues(dofLocalNo, values);
}

//! for a specific component, get a single value from local dof no.
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofLocalNo)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  return this->component_[componentNo].getValue(dofLocalNo);
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nElementsLocal());
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getElementValues(elementNo, values);
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nElementsLocal());
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();

  // TODO: local to global
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
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  this->values_ = rhs.partitionedPetscVec();
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->nComponents == 1 && nComponents == 1);
  const int nValues = values.size();

  this->values_->setValues(0, nValues, (const int *) dofNosLocal.data(), values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set value for all dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(double value)
{
  // get number of dofs
  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nDofsLocal();

  std::vector<double> valueBuffer(nDofs,value);
  
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->setValues(componentIndex, valueBuffer, INSERT_VALUES);
  }
}

//! set values for all components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(dofNosLocal.size() == values.size());
 
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

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set a single dof (all components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode)
{
  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->setValues(componentIndex, 1, &dofLocalNo, value.data()+componentIndex, petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValuesWithGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithGhosts());
 
  this->values_->setValues(componentNo, values.size(), this->mesh_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValuesWithoutGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithoutGhosts());
 
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->values_->setValues(componentNo, values.size(), this->mesh_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

//! set values for all components for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValuesWithGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithGhosts());
  
  this->setValues(this->mesh_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

//! set values for all components for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValuesWithoutGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithoutGhosts());
  
  this->setValues(this->mesh_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

//! set value to zero for all dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
zeroEntries()
{
  this->values_->zeroEntries();
}
  
//! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
finishVectorManipulation()
{
  this->values_->finishVectorManipulation();
}

};  // namespace
