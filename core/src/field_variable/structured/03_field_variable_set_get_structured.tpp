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
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  assert(componentNo >= 0 && componentNo < nComponents);
 
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->mesh_->meshPartition()->nDofsLocalWithGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }
  
  // resize output vector
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << nValues << " entries.";
  values.resize(nValues);
  
  // get values
  this->values_->getValues(componentNo, nValues, this->mesh_->meshPartition()->dofNosLocal(onlyNodalValues).data(), values.data());
}

//! for a specific component, get all values
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  assert(componentNo >= 0 && componentNo < nComponents);
 
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->mesh_->meshPartition()->nDofsLocalWithGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }
  
  // resize output vector
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << previousSize+nValues << " entries.";
  
  // get values
  this->values_->getValues(componentNo, nValues, this->mesh_->meshPartition()->dofNosLocal(onlyNodalValues).data(), values.data()+previousSize);
}

//! for a specific component, get values from their local dof no.s
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofLocalNo, std::vector<double> &values)
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
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->values_->getValues(componentNo, N, (PetscInt *)dofLocalNo.data(), values.data());
}

template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values)
{
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

//! for a specific component, get the values corresponding to all element-local dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo,
                 std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nElementsLocal());
  assert(componentNo >= 0 && componentNo < nComponents);
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();

  // create indices with dofs
  std::array<PetscInt,nDofsPerElement> indices;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = this->mesh_->getDofNo(elementNo, dofIndex);
  }
  
  // get the values
  this->values_->getValues(componentNo, nDofsPerElement, indices.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, 
                 std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nElementsLocal());
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  std::array<PetscInt,nDofsPerElement> indices;
  std::array<double,nDofsPerElement*nComponents> result;

  VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents << ", nDofsPerElement=" << nDofsPerElement;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      indices[dofIndex] = this->mesh_->getDofNo(elementNo, dofIndex);
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
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofLocalNo)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  PetscInt index = dofLocalNo;

  double result;
  this->values_->getValues(componentNo, 1, &index, &result);
  return result;
}
  
//! set values for all components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  if (dofNosLocal.size() != values.size())
  {
    LOG(DEBUG) << "dofNosLocal.size(): " << dofNosLocal.size() << ", values.size(): " << values.size();
    LOG(DEBUG) << "dofNosLocal: " << dofNosLocal << ", values: " << values;
  }
  assert(dofNosLocal.size() == values.size());

  setValues(values.size(), dofNosLocal, values, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for all components for dofs, only nValues values will be set despite potentially more dofNosLocal, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValues(int nValues, const std::vector<dof_no_t> &dofNosLocal, const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(dofNosLocal.size() >= nValues);
  assert(values.size() == nValues);

  std::vector<double> valuesBuffer(nValues);

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

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set a single dof (all components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode)
{
  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->setValues(componentIndex, 1, &dofLocalNo, value.data()+componentIndex, petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set value for all dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValues(double value)
{
  // get number of dofs
  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nDofsLocalWithGhosts();

  std::vector<double> valueBuffer(nDofs,value);
  
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->setValuesWithGhosts(componentIndex, valueBuffer, INSERT_VALUES);
  }
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValuesWithGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithGhosts());
 
  // set the values
  this->values_->setValues(componentNo, values.size(), this->mesh_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

//! set values for the specified component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValuesWithoutGhosts(int componentNo, const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithoutGhosts());
 
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->values_->setValues(componentNo, values.size(), this->mesh_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValuesWithGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithGhosts());
  
  this->setValues(this->mesh_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValuesWithoutGhosts(const std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithoutGhosts());
  
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->setValues(values.size(), this->mesh_->meshPartition()->dofNosLocal(), values, petscInsertMode);
}

//! set value to zero for all dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
zeroEntries()
{
  this->values_->zeroEntries();
}
  
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
finishVectorManipulation()
{
  this->values_->finishVectorManipulation();
}

};
