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
getValues(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  assert(componentNo >= 0 && componentNo < nComponents);
 
  // set stride to nDofsPerNode if Hermite, else to 1
  const int stride = (onlyNodalValues && std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value ? BasisOnMeshType::nDofsPerNode() : 1);
  
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = nDofs;
  if (onlyNodalValues)
    // if the basis function is Hermite
    if (std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value)
      nValues = nDofs / BasisOnMeshType::nDofsPerNode();

  // store the array indices for values_ array in dofLocalNo
  std::vector<PetscInt> indices(nValues,0);
  dof_no_t indexNo = 0;
  for (dof_no_t dofLocalNo = 0; dofLocalNo < nDofs; dofLocalNo += stride)
  {
    assert(indexNo < nValues);
    indices[indexNo++] = dofLocalNo;
  }

  VLOG(2) << "Field variable structured, getValues, resize values vector to " << nValues << " entries.";
  values.resize(nValues);
  this->values_->getValues(componentNo, nValues, indices.data(), values.data());
}

//! for a specific component, get values from their local dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
 
  // store the array indices for values_ array in dofLocalNo
  for (int i=0; i<N; i++)
  {
    dofLocalNo[i] = dofLocalNo[i];
  }

  this->values_->getValues(componentNo, N, (PetscInt *)dofLocalNo.data(), values.data());
}

//! for a specific component, get values from their local dof no.s
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofLocalNo, std::vector<double> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  // store the array indices for values_ array in dofLocalNo
  const int nValues = dofLocalNo.size();
  for (int i=0; i<nValues; i++)
  {
    dofLocalNo[i] = dofLocalNo[i];
    
    assert(dofLocalNo[i] >= 0 && dofLocalNo[i] < this->nEntries_);
  }

  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VLOG(3) << "FieldVariable::getValues VecGetValues(Vec@" << this->values_ << "," << nValues << "," << (PetscInt *)dofLocalNo.data() << "," << values.data()+previousSize << "," << previousSize << ")";
  
  this->values_->getValues(componentNo, nValues, (PetscInt *)dofLocalNo.data(), values.data()+previousSize);
}

template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values)
{
  std::array<double,N*nComponents> result;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
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
  assert(elementNo >= 0 && elementNo < this->mesh_->nLocalElements());
  assert(componentNo >= 0 && componentNo < nComponents);
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();

  std::array<PetscInt,nDofsPerElement> indices;

  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = this->mesh_->getDofNo(elementNo,dofIndex);
  }

  this->values_->getValues(componentNo, nDofsPerElement, indices.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nLocalElements());
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  std::array<int,nComponents> indices;
  std::array<double,nDofsPerElement*nComponents> result;

  //VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
  {
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      indices[dofIndex] = this->mesh_->getDofNo(elementNo,dofIndex);
    }
    
    this->values_->getValues(componentIndex, nDofsPerElement, indices.data(), result.data() + componentIndex*nDofsPerElement);
  }

  //VLOG(2) << " indices: " << indices << ", retrieved values: " << result;

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
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
setValues(std::vector<dof_no_t> &dofLocalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  assert(dofLocalNos.size() == values.size());
 
  const int nValues = values.size();
  std::array<double,nValues> valuesBuffer;

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // loop over dofs and prepare values of current component
    for (int dofIndex = 0; dofIndex < nValues; dofIndex++)
    {
      valuesBuffer[dofIndex] = values[dofIndex][componentIndex];
    }
    
    this->values_->setValues(componentIndex, nValue, dofLocalNos.data(), valuesBuffer.data(), petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set a single dof (all components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValue(dof_no_t dofLocalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode)
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
  const dof_no_t nDofs = this->mesh_->nLocalDofs();

  std::array<PetscInt, nDofs> indices;
  std::array<double, nDofs> valueBuffer;
  valueBuffer.fill(value);
  
  std::iota(indices.begin(), indices.end(), 0);
  
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    this->values_->setValues(componentIndex, nDofs, indices.data(), valueBuffer.data(), INSERT_VALUES);
  }
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
finishVectorManipulation()
{
  this->values_->finishVectorManipulation();
}

};
