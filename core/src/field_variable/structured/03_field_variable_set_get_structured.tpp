#include "field_variable/field_variable.h"

#include <sstream>
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
  
  const dof_no_t nDofs = this->mesh_->nDofs();

  // set stride to nDofsPerNode if Hermite, else to 1
  const int stride = (onlyNodalValues && std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value ? BasisOnMeshType::nDofsPerNode() : 1);
  
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = nDofs;
  if (onlyNodalValues)
    // if the basis function is Hermite
    if (std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value)
      nValues = nDofs / BasisOnMeshType::nDofsPerNode();
    
  VLOG(2) << "getValues";
  VLOG(2) << "componentNumber=" << componentNo << ", nDofs=" << nDofs << ", nValues="<< nValues << ", stride=" << stride << ", nComponents=" << nComponents;
  
  // store the array indices for values_ array in dofGlobalNo
  std::vector<PetscInt> indices(nValues,0);
  dof_no_t indexNo = 0;
  for (dof_no_t dofGlobalNo=0; dofGlobalNo<nDofs; dofGlobalNo+=stride)
  {
    assert(indexNo < nValues);
    indices[indexNo++] = componentNo*nDofs + dofGlobalNo;
  
    VLOG(2) << "Indices of " << indexNo-1 << ": " << indices[indexNo-1];
    
  }

  VLOG(2) << "Field variable structured, getValues, resize values vector to " << nValues << " entries.";
  values.resize(nValues);
  VecGetValues(this->values_, nValues, indices.data(), values.data());
  
  if (VLOG_IS_ON(2))
  {
    VLOG(2) << "Retrieved values: " <<  values;
    VLOG(2) << "Petsc vector: " << PetscUtility::getStringVector(this->values_);
  }
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  const dof_no_t nDofs = this->mesh_->nDofs();
  
  // store the array indices for values_ array in dofGlobalNo
  for (int i=0; i<N; i++)
  {
    dofGlobalNo[i] = componentNo*nDofs + dofGlobalNo[i];
  }

  VecGetValues(this->values_, N, (PetscInt *)dofGlobalNo.data(), values.data());
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  const dof_no_t nDofs = this->mesh_->nDofs();
  
  // store the array indices for values_ array in dofGlobalNo
  const int nValues = dofGlobalNo.size();
  for (int i=0; i<nValues; i++)
  {
    dofGlobalNo[i] = componentNo*nDofs + dofGlobalNo[i];
    
    assert(dofGlobalNo[i] >= 0 && dofGlobalNo[i] < this->nEntries_);
  }

  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VLOG(3) << "FieldVariable::getValues VecGetValues(Vec@" << this->values_ << "," << nValues << "," << (PetscInt *)dofGlobalNo.data() << "," << values.data()+previousSize << "," << previousSize << ")";
  
  VecGetValues(this->values_, nValues, (PetscInt *)dofGlobalNo.data(), values.data()+previousSize);
}

template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  std::array<int,N*nComponents> indices;
  std::array<double,N*nComponents> result;

  const dof_no_t nDofs = this->mesh_->nDofs();
  
  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int i=0; i<N; i++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      indices[j] = componentIndex*nDofs + dofGlobalNo[i];
    }
  }

  VecGetValues(this->values_, N*nComponents, indices.data(), result.data());

  // copy result to output values
  for (int i=0; i<N; i++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      values[i][componentIndex] = result[i*nComponents+componentIndex];
    }
  }
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo,
                 std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nElements());
  assert(componentNo >= 0 && componentNo < nComponents);
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const dof_no_t nDofs = this->mesh_->nDofs();
  
  std::array<PetscInt,nDofsPerElement> indices;

  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = componentNo*nDofs + this->mesh_->getDofNo(elementNo,dofIndex);
  }

  VecGetValues(this->values_, nDofsPerElement, indices.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nElements());
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const dof_no_t nDofs = this->mesh_->nDofs();
  
  std::array<int,nDofsPerElement*nComponents> indices;
  std::array<double,nDofsPerElement*nComponents> result;

  VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents;

  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      indices[j] = componentIndex*nDofs + this->mesh_->getDofNo(elementNo,dofIndex);
    }
  }

  VecGetValues(this->values_, nDofsPerElement*nComponents, indices.data(), result.data());

  if (VLOG_IS_ON(2))
  {
    VLOG(2) << "Indices: " << indices << ", retrieved values: " << result;
    VLOG(2) << "Values: " << PetscUtility::getStringVector(this->values_);
  }

  // copy result to output values
  for (int dofIndex=0; dofIndex<nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      values[dofIndex][componentIndex] = result[dofIndex*nComponents+componentIndex];
      VLOG(2) << "getElementValues element " << elementNo << ", dofIndex " << dofIndex << " componentIndex " << componentIndex << " value: " << values[dofIndex][componentIndex];
    }
  }
}

//! for a specific component, get a single value from global dof no.
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofGlobalNo)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  const dof_no_t nDofs = this->mesh_->nDofs();
  
  PetscInt index = componentNo*nDofs + dofGlobalNo;

  double result;
  VecGetValues(this->values_, 1, &index, &result);
  return result;
}

//! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  std::array<int,nComponents> indices;

  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nDofs();
  
  // loop over dof numbers
  int i=0;
  for (std::vector<dof_no_t>::iterator iter = dofGlobalNos.begin(); iter != dofGlobalNos.end(); iter++, i++)
  {
    dof_no_t dofGlobalNo = *iter;

    // prepare lookup indices for PETSc vector values_
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = componentIndex*nDofs + dofGlobalNo;
    }

    VecSetValues(this->values_, nComponents, indices.data(), values[i].data(), petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

//! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode)
{
  std::array<int,nComponents> indices;
  const dof_no_t nDofs = this->mesh_->nDofs();

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    indices[componentIndex] = componentIndex*nDofs + dofGlobalNo;
  }

  VecSetValues(this->values_, nComponents, indices.data(), value.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

//! set value for all dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
setValues(double value)
{
  VecSet(this->values_, value);
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
flushSetValues()
{
  VecAssemblyBegin(this->values_);
  VecAssemblyEnd(this->values_);
}

};
