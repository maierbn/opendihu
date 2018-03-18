#include "field_variable/field_variable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>
#include <array>

namespace FieldVariable
{
//! for a specific component, get all values
template<typename BasisOnMeshType>
void FieldVariableSetGetStructured<BasisOnMeshType>::
getValues(std::string component, std::vector<double> &values)
{
  int componentIndex = this->componentIndex_[component];
  const dof_no_t nDofs = this->mesh_->nDofs();
  values.resize(nDofs);
  
  // store the array indices for values_ array in dofGlobalNo
  std::vector<int> indices(nDofs,0);
  for (dof_no_t dofGlobalNo=0; dofGlobalNo<nDofs; dofGlobalNo++)
  {
    indices[dofGlobalNo] = dofGlobalNo*this->nComponents_ + componentIndex;
  }
  
  VecGetValues(this->values_, nDofs, indices.data(), values.data());
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType>::
getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  int componentIndex = this->componentIndex_[component];
 
  // store the array indices for values_ array in dofGlobalNo
  for (int i=0; i<N; i++)
  {
    dofGlobalNo[i] = dofGlobalNo[i]*this->nComponents_ + componentIndex;
  }
  
  VecGetValues(this->values_, N, dofGlobalNo.data(), values.data());
}

/*
//! get values from their global dof no.s for all components
template<typename BasisOnMeshType>
template<int N, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  std::array<int,N*nComponents> indices;
  std::array<double,N*nComponents> result;
  
  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int i=0; i<N; i++)
  {
    for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
    {
      indices[j] = dofGlobalNo[i]*nComponents + componentIndex;
    }
  }
  
  VecGetValues(this->values_, N*nComponents, indices.data(), result.data());
  
  // copy result to output values
  for (int i=0; i<N; i++)
  {
    for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
    {
      values[i][componentIndex] = result[i*nComponents+componentIndex];
    }
  }
}
*/

//! for a specific component, get the values corresponding to all element-local dofs
template<typename BasisOnMeshType>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType>::
getElementValues(std::string component, element_no_t elementNo, 
                 std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  int componentIndex = this->componentIndex_[component];
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  std::array<int,nDofsPerElement> indices;
  
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = this->mesh_->getDofNo(elementNo,dofIndex)*this->nComponents_ + componentIndex;
  }
  
  VecGetValues(this->values_, N, indices.data(), values.data());
}
/*

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType>
template<std::size_t nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  std::array<int,nDofsPerElement*nComponents> indices;
  std::array<double,nDofsPerElement*nComponents> result;
  
  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
    {
      indices[j] = this->mesh_->getDofNo(elementNo,dofIndex)*nComponents + componentIndex;
    }
  }
  
  VecGetValues(this->values_, nDofsPerElement*nComponents, indices.data(), result.data());
  
  // copy result to output values
  for (int dofIndex=0; dofIndex<nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
    {
      values[dofIndex][componentIndex] = result[dofIndex*nComponents+componentIndex];
    }
  }
}
*/

//! for a specific component, get a single value from global dof no.
template<typename BasisOnMeshType>
double FieldVariableSetGetStructured<BasisOnMeshType>::
getValue(std::string component, node_no_t dofGlobalNo)
{
  int componentIndex = this->componentIndex_[component];
  int index = dofGlobalNo*this->nComponents_ + componentIndex;
  
  double result;
  VecGetValues(this->values_, 1, &index, &result); 
  return result;
}

//! get a single value from global dof no. for all components
template<typename BasisOnMeshType>
template<std::size_t nComponents>
std::array<double,nComponents> FieldVariableSetGetStructured<BasisOnMeshType>::
getValue(node_no_t dofGlobalNo)
{
  std::array<int,nComponents> indices;
  std::array<double,nComponents> result;
  
  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
  {
    indices[j] = dofGlobalNo*this->nComponents_ + componentIndex;
  }
    
  VecGetValues(this->values_, nComponents, indices.data(), result.data());
  return result;
}

/*
//! set values for dofs
template<typename BasisOnMeshType>
template<std::size_t nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values)
{
  std::array<int,nComponents> indices;
  
  // loop over dof numbers
  int i=0;
  for (std::vector<dof_no_t>::iterator iter = dofGlobalNos.begin(); iter != dofGlobalNos.end(); iter++, i++)
  {  
    dof_no_t dofGlobalNo = *iter;
    
    // prepare lookup indices for PETSc vector values_
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = dofGlobalNo*this->nComponents_ + componentIndex;
    }
    
    VecSetValues(this->values_, nComponents, indices.data(), values[i].data(), INSERT_VALUES);
  }
  
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}*/

/*
//! set a single value
template<typename BasisOnMeshType>
template<std::size_t nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType>::
setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value)
{
  std::array<int,nComponents> indices;
  
  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    indices[componentIndex] = dofGlobalNo*this->nComponents_ + componentIndex;
  }
  
  VecSetValues(this->values_, nComponents, indices.data(), value.data(), INSERT_VALUES);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}*/

template<typename BasisOnMeshType>
void FieldVariableSetGetStructured<BasisOnMeshType>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->nComponents == 1);
  const int nValues = values.size();

  VecSetValues(this->values_, nValues, (const int *) dofGlobalNos.data(), values.data(), petscInsertMode);
  
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}

template<typename BasisOnMeshType>
void FieldVariableSetGetStructured<BasisOnMeshType>::
flushSetValues()
{
  VecAssemblyBegin(this->values_); 
  VecAssemblyEnd(this->values_);
}

};
