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
  assert(this->values_);

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
  assert(this->values_);

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

//! get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValuesWithGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues) const
{
  assert(this->values_);

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

  // loop over components and get data component-wise
  std::vector<double> buffer(nValues);
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // get values into buffer
    this->values_->getValues(componentNo, nValues, this->functionSpace_->meshPartition()->dofNosLocal(onlyNodalValues).data(), buffer.data());

    // copy values from buffer to output vector
    for (int valueIndex = 0; valueIndex < nValues; valueIndex++)
    {
      values[valueIndex][componentNo] = buffer[valueIndex];
    }
  }
}

//! get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValuesWithoutGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues) const
{
  assert(this->values_);

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

  // loop over components and get data component-wise
  std::vector<double> buffer(nValues);
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // get values into buffer
    this->values_->getValues(componentNo, nValues, this->functionSpace_->meshPartition()->dofNosLocal(onlyNodalValues).data(), buffer.data());

    // copy values from buffer to output vector
    for (int valueIndex = 0; valueIndex < nValues; valueIndex++)
    {
      values[previousSize+valueIndex][componentNo] = buffer[valueIndex];
    }
  }
}

//! get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValuesWithoutGhosts(std::array<std::vector<double>,nComponents> &values, bool onlyNodalValues) const
{
  assert(this->values_);

  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }

  // loop over components and get data component-wise
  std::vector<double> buffer(nValues);
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::size_t previousSize = values[componentNo].size();
    values[componentNo].resize(previousSize+nValues);

    // resize output vector
    VLOG(2) << "Field variable structured, getValues, resize values vector to " << previousSize+nValues << " entries.";

    // get values into buffer
    this->values_->getValues(componentNo, nValues, this->functionSpace_->meshPartition()->dofNosLocal(onlyNodalValues).data(), values[componentNo].data());
  }
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
getValues(int componentNo, int nValues, const dof_no_t *dofLocalNo, std::vector<double> &values) const
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  int valuesPreviousSize = values.size();
  values.resize(valuesPreviousSize + nValues);
  this->values_->getValues(componentNo, nValues, (PetscInt *)dofLocalNo, values.data() + valuesPreviousSize);
}

//! get values from their local dof no.s for all components
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValues(std::vector<dof_no_t> dofLocalNo, std::vector<std::array<double,nComponents>> &values) const
{
  assert(this->values_);
  const int nValues = dofLocalNo.size();
  std::vector<double> result(nValues*nComponents);   // temporary result buffer

  int initialSize = values.size();
  values.resize(initialSize + nValues);

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->getValues(componentIndex, nValues, dofLocalNo.data(), result.data() + componentIndex*nValues);
  }

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nValues; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      values[initialSize+dofIndex][componentIndex] = result[componentIndex*nValues + dofIndex];
    }
  }
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
{
  assert(this->values_);
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
getElementValues(int componentNo, element_no_t elementNoLocal,
                 std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(elementNoLocal >= 0 && elementNoLocal < this->functionSpace_->nElementsLocal());
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);
  
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // create indices with dofs
  std::array<PetscInt,nDofsPerElement> indices;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = this->functionSpace_->getDofNo(elementNoLocal, dofIndex);
    //LOG(DEBUG) << "getElementValues el. " << elementNoLocal << ", dof " << indices[dofIndex];
  }
  
  // get the values
  this->values_->getValues(componentNo, nDofsPerElement, indices.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getElementValues(element_no_t elementNoLocal,
                 std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(this->functionSpace_);
#ifndef NDEBUG
  if (elementNoLocal >= this->functionSpace_->nElementsLocal())
  {
    LOG(ERROR) << "getElementValues, elementNoLocal = " << elementNoLocal << " >= nElementsLocal = " << this->functionSpace_->nElementsLocal();
  }
#endif
  assert(elementNoLocal >= 0 && elementNoLocal < this->functionSpace_->nElementsLocal());
  assert(this->values_);

  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  std::array<PetscInt,nDofsPerElement> indices;
  std::array<double,nDofsPerElement*nComponents> result;

  VLOG(2) << "getElementValues element " << elementNoLocal << ", nComponents=" << nComponents << ", nDofsPerElement=" << nDofsPerElement;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      indices[dofIndex] = this->functionSpace_->getDofNo(elementNoLocal, dofIndex);
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
      //VLOG(2) << "getElementValues element " << elementNoLocal << ", dofIndex " << dofIndex << " componentIndex " << componentIndex << " value: " << values[dofIndex][componentIndex];
    }
  }
}

//! vectorized version of getElementValues
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getElementValues(Vc::int_v elementNoLocal,
                 std::array<std::array<Vc::double_v,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(this->functionSpace_);
  assert(this->values_);

  const int nVcComponents = Vc::double_v::size();
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  std::array<PetscInt, nDofsPerElement*nVcComponents> indices;
  std::array<Vc::double_v, nDofsPerElement*nComponents> result;

  VLOG(2) << "getElementValues (vectorized) element " << elementNoLocal << ", nComponents=" << nComponents << ", nDofsPerElement=" << nDofsPerElement << ", nVcComponents: " << nVcComponents;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int vcComponent = 0; vcComponent < nVcComponents; vcComponent++)
      {
        if (elementNoLocal[vcComponent] == -1)
          indices[dofIndex*nVcComponents + vcComponent] = 0;    // set index to 0, then here the dof 0 is retrieved and also further used in computation but it is discarded later in setValue
        else
          indices[dofIndex*nVcComponents + vcComponent] = this->functionSpace_->getDofNo(elementNoLocal[vcComponent], dofIndex);
        //LOG(DEBUG) << "    element " << elementNoLocal[vcComponent] << " dof " << dofIndex << ": indices[" << dofIndex*nVcComponents + vcComponent << "]: " << indices[dofIndex*nVcComponents + vcComponent];
      }
    }

    // get the values for the current component
    this->values_->getValues(componentIndex, nDofsPerElement*nVcComponents, indices.data(), (double *)&result[componentIndex*nDofsPerElement]);

    //LOG(DEBUG) << "component " << componentIndex << ", get " << nDofsPerElement << " values at " << indices << ", result (starting at index " << componentIndex*nDofsPerElement << "): " << result;

  }

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      values[dofIndex][componentIndex] = result[componentIndex*nDofsPerElement + dofIndex];

      //LOG(DEBUG) << "dof " << dofIndex << " component " << componentIndex << " value " << result[componentIndex*nDofsPerElement + dofIndex] << " -> " << values[dofIndex][componentIndex];
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

//! get a single value from local dof no. for all components
template<typename FunctionSpaceType, int nComponents>
std::array<double,nComponents> FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
getValue(node_no_t dofLocalNo) const
{
  assert(this->values_);

  PetscInt index = dofLocalNo;

  std::array<double,nComponents> result;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    this->values_->getValues(componentNo, 1, &index, &result[componentNo]);
  }
  return result;
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
extractComponentCopy(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable)
{
  assert(extractedFieldVariable->partitionedPetscVec());
  assert(this->values_);
  this->values_->extractComponentCopy(componentNo, extractedFieldVariable->partitionedPetscVec());
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
extractComponentShared(int componentNo, std::shared_ptr<FieldVariable<FunctionSpaceType,1>> extractedFieldVariable)
{
  assert(extractedFieldVariable->partitionedPetscVec());
  assert(this->values_);
  this->values_->extractComponentShared(componentNo, extractedFieldVariable->partitionedPetscVec());
}

template<typename FunctionSpaceType, int nComponents>
bool FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
isExtractComponentSharedPossible(int componentNo)
{
  assert(this->values_);
  return this->values_->isExtractComponentSharedPossible(componentNo);
}

template<typename FunctionSpaceType, int nComponents>
template<int nComponents2>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
restoreExtractedComponent(std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents2>> extractedVec, int componentNo)
{
  assert(this->values_);
  this->values_->template restoreExtractedComponent<nComponents2>(extractedVec, componentNo);
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

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int componentNo, const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);
  if (values.size() < dofNosLocal.size())
  {
    LOG(FATAL) << "FieldVariable::setValues: trying to set " << dofNosLocal.size() << " values but only " << values.size() << " given.";
  }
  assert(values.size() >= dofNosLocal.size());

  // set the values for the given component
  this->values_->setValues(componentNo, dofNosLocal.size(), dofNosLocal.data(), values.data(), petscInsertMode);
}

template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int componentNo, const std::array<dof_no_t,N> &dofNosLocal, const std::array<double,N> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  // set the values for the given component
  this->values_->setValues(componentNo, N, dofNosLocal.data(), values.data(), petscInsertMode);
}

//! set values for a given component for given dofs, using raw pointers
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  // set the values for the given component
  this->values_->setValues(componentNo, nValues, dofNosLocal, values, petscInsertMode);
}

//! set values for all components for dofs, only nValues values will be set despite potentially more dofNosLocal, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(int nValues, const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  if (nValues == 0)
    return;

  assert(this->values_);
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

  // loop over components and set single value for each component
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->setValues(componentIndex, 1, &dofLocalNo, value.data()+componentIndex, petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set a single dof (all components), after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValue(Vc::int_v dofLocalNo, const std::array<Vc::double_v,nComponents> &value, InsertMode petscInsertMode)
{
  // loop over components and set vectorized value for each component
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->setValue(componentIndex, dofLocalNo, value[componentIndex], petscInsertMode);
  }
}

//! set a given component of Vc::double_v::size() dofs with the vectorized value, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValue(int componentNo, Vc::int_v dofLocalNo, Vc::double_v value, InsertMode petscInsertMode)
{
  assert(this->values_);

  // count number of non-negative indices in dofLocalNo, it is assumed that they occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::isnegative(dofLocalNo).count();
/*
  // store Vc vectors in order to get the raw memory
  std::array<double,Vc::double_v::size()> data;
//  for (int i = 0; i < Vc::double_v::size(); i++)
//    data[i] = value[i];
  
  value.store(data.data(),Vc::Aligned);

  std::array<int,Vc::int_v::size()> indices;
  dofLocalNo.store(indices.data());
  VLOG(1) << "setValue(componentNo=" << componentNo << ", dofLocalNo=" << dofLocalNo << ", value=" << value << ")";
  VLOG(1) << "indices: " << indices << ", data: " << data;
*/
  this->values_->setValues(componentNo, nEntries, (PetscInt *)&dofLocalNo, (double *)&value, petscInsertMode);
}

//! set a given component of Vc::double_v::size() dofs with the same value
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValue(int componentNo, Vc::int_v dofLocalNo, double value, InsertMode petscInsertMode)
{
  assert(this->values_);
  std::array<double,Vc::double_v::size()> data;
  data.fill(value);

  // store Vc vectors in order to get the raw memory
  std::array<int,Vc::int_v::size()> indices;
  dofLocalNo.store(indices.data());

  // count number of non-negative indices in dofLocalNo, it is assumed that they occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::isnegative(dofLocalNo).count();

  this->values_->setValues(componentNo, nEntries, indices.data(), data.data(), petscInsertMode);
}

template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValues(const std::array<dof_no_t,N> &dofNosLocal, const std::array<std::array<double,nComponents>,N> &values, InsertMode petscInsertMode)
{
  assert(this->values_);

  std::array<double,N> valuesBuffer;

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // loop over dofs and prepare values of current component
    for (int dofIndex = 0; dofIndex < N; dofIndex++)
    {
      valuesBuffer[dofIndex] = values[dofIndex][componentIndex];
    }

    // set the values for the current component
    this->values_->setValues(componentIndex, N, dofNosLocal.data(), valuesBuffer.data(), petscInsertMode);
  }
}

//! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValue(int componentNo, dof_no_t dofLocalNo, double value, InsertMode petscInsertMode)
{
  assert(this->values_);

  this->values_->setValues(componentNo, 1, &dofLocalNo, &value, petscInsertMode);
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

template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
setValuesWithoutGhosts(const std::array<std::vector<double>,nComponents> &values, InsertMode petscInsertMode)
{
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    if (values[componentIndex].size() != this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts())
    {
      LOG(FATAL) << "setValuesWithoutGhosts: values.size: " << values[componentIndex].size() << " != " << this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
    }
    assert(values[componentIndex].size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());

    // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
    this->setValues(componentIndex, this->functionSpace_->meshPartition()->dofNosLocal(), values[componentIndex], petscInsertMode);
  }
}

//! set value to zero for all dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
zeroEntries()
{
  assert(this->values_);
  this->values_->zeroEntries();
}
}  // namespace
