#include "field_variable/field_variable.h"

#include <sstream>
#include <algorithm>
#include <cassert>
#include <array>
#include <type_traits>

#include "basis_function/hermite.h"

namespace FieldVariable
{
//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
 
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->meshPartition()->nDofsLocalWithGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }
  
  // resize output vector
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << nValues << " entries.";
  values.resize(nValues);
  
  // get values
  this->values_->getValues(componentNo, nValues, this->functionSpace_->meshPartition()->dofNosLocal(onlyNodalValues).data(), values.data());
}

//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
 
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }
  
  // resize output vector
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << previousSize+nValues << " entries.";
  
  // get values
  this->values_->getValues(componentNo, nValues, this->functionSpace_->meshPartition()->dofNosLocal(onlyNodalValues).data(), values.data()+previousSize);
}

//! for a specific component, get values from their local dof no.s
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValues(int componentNo, const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  // store the array indices for values_ array in dofLocalNo
  const int nValues = dofLocalNo.size();

  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << previousSize+nValues << " entries.";
  
  this->values_->getValues(componentNo, nValues, (PetscInt *)dofLocalNo.data(), values.data()+previousSize);
}

//! for a specific component, get values from their local dof no.s
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  this->values_->getValues(componentNo, N, (PetscInt *)dofLocalNo.data(), values.data());
}

template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values) const
{
  assert(this->values_);
  std::array<double,N*nComponents> result;   // temporary result buffer

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->getValues(componentIndex, N, dofLocalNo.data(), result.data() + componentIndex*N);
  }

  // copy result to output values
  for (int dofIndex = 0; dofIndex < N; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      values[dofIndex][componentIndex] = result[componentIndex*N + dofIndex];
    }
  }
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
{
  int nValues = dofLocalNo.size();
  values.resize(nValues*nComponents);

  // get values into values vector
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->getValues(componentIndex, nValues, dofLocalNo.data(), values.data() + componentIndex*nValues);
  }
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo,
                 std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);
  
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // create indices with dofs
  std::array<PetscInt,nDofsPerElement> indices;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = this->functionSpace_->getDofNo(elementNo, dofIndex);
  }
  
  // get the values
  this->values_->getValues(componentNo, nDofsPerElement, indices.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getElementValues(element_no_t elementNo, 
                 std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(this->values_);
  
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  std::array<PetscInt,nDofsPerElement> indices;
  std::array<double,nDofsPerElement*nComponents> result;

  VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents << ", nDofsPerElement=" << nDofsPerElement;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      indices[dofIndex] = this->functionSpace_->getDofNo(elementNo, dofIndex);
    }
    
    // get the values for the current component
    this->values_->getValues(componentIndex, nDofsPerElement, indices.data(), result.data() + componentIndex*nDofsPerElement);
  }

  //VLOG(2) << " indices: " << indices << ", retrieved values: " << result;

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      values[dofIndex][componentIndex] = result[componentIndex*nDofsPerElement + dofIndex];
      //VLOG(2) << "getElementValues element " << elementNo << ", dofIndex " << dofIndex << " componentIndex " << componentIndex << " value: " << values[dofIndex][componentIndex];
    }
  }
}

//! for a specific component, get a single value from local dof no.
template<typename FunctionSpaceType, int nComponents>
double FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValue(int componentNo, node_no_t dofLocalNo) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);
  
  PetscInt index = dofLocalNo;

  double result;
  this->values_->getValues(componentNo, 1, &index, &result);
  return result;
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
extractComponentCopy(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable)
{
  assert(extractedFieldVariable->partitionedPetscVec());
  this->values_->extractComponentCopy(componentNo, extractedFieldVariable->partitionedPetscVec());
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
extractComponentShared(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable)
{
  assert(extractedFieldVariable->partitionedPetscVec());
  this->values_->extractComponentShared(componentNo, extractedFieldVariable->partitionedPetscVec());
}

template<typename FunctionSpaceType, int nComponents>
template<int nComponents2>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
restoreExtractedComponent(std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents2>> extractedVec)
{
  this->values_->template restoreExtractedComponent<nComponents2>(extractedVec);
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int componentNo, Vec petscVector)
{
  this->values_->setValues(componentNo, petscVector);
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> fieldVariable)
{
  assert(fieldVariable->partitionedPetscVec());
  this->values_->setValues(componentNo, fieldVariable->partitionedPetscVec());
}

//! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
#ifndef NDEBUG
  if (dofNosLocal.size() != values.size())
  {
    LOG(DEBUG) << "dofNosLocal.size(): " << dofNosLocal.size() << ", values.size(): " << values.size();
    LOG(DEBUG) << "dofNosLocal: " << dofNosLocal << ", values: " << values;
  }
  assert(dofNosLocal.size() == values.size());
#endif

  setValues(values.size(), dofNosLocal, values, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set values for all components for dofs, only nValues values will be set despite potentially more dofNosLocal, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int nValues, const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  if (nValues == 0)
    return;

  if (dofNosLocal.size() < nValues)
  {
    LOG(ERROR) << "dofNosLocal.size()=" << dofNosLocal.size() << ", nValues=" << nValues;
  }
  assert(dofNosLocal.size() >= nValues);
  assert(values.size() == nValues);
  assert(this->values_);

  static std::vector<double> valuesBuffer;
  valuesBuffer.resize(nValues);

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // loop over dofs and prepare values of current component
    for (int dofIndex = 0; dofIndex < nValues; dofIndex++)
    {
      valuesBuffer[dofIndex] = values[dofIndex][componentIndex];
    }

    // set the values for the current component
    this->values_->setValues(componentIndex, nValues, dofNosLocal.data(), valuesBuffer.data(), petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
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

//! set value for all dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(double value)
{
  // get number of dofs
  assert(this->functionSpace_);
  const dof_no_t nDofs = this->functionSpace_->nDofsLocalWithGhosts();

  static std::vector<double> valueBuffer;
  valueBuffer.assign(nDofs,value);
  
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->setValuesWithGhosts(componentIndex, valueBuffer, INSERT_VALUES);
  }
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValuesWithGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());
  assert(this->values_);
 
  // set the values
  this->values_->setValues(componentNo, values.size(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValuesWithoutGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
  assert(this->values_);
 
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->values_->setValues(componentNo, values.size(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValuesWithGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());
  
  this->setValues(this->functionSpace_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValuesWithoutGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  if (values.size() != this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts())
  {
    LOG(FATAL) << "setValuesWithoutGhosts: values.size: " << values.size() << " != " << this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  }
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
  
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->setValues(values.size(), this->functionSpace_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

//! set value to zero for all dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
zeroEntries()
{
  assert(this->values_);
  this->values_->zeroEntries();
}
};
