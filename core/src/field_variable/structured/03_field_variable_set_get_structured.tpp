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
  const dof_no_t nDofs = this->mesh_->nDofs();

  // set stride to 2 if Hermite, else to 1  
  const int stride = (onlyNodalValues && std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value ? 2 : 1);
  
  // determine the number of values to be retrived which is half the number of dofs for Hermite with only nodal values
  dof_no_t nValues = nDofs;
  if (onlyNodalValues)
    // if the basis function is Hermite
    if (std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value)
      nValues = nDofs / 2;
  
  // store the array indices for values_ array in dofGlobalNo
  std::vector<PetscInt> indices(nValues,0);
  dof_no_t indexNo = 0;
  for (dof_no_t dofGlobalNo=0; dofGlobalNo<nDofs; dofGlobalNo+=stride)
  {
    assert(indexNo < nValues);
    indices[indexNo++] = dofGlobalNo*nComponents + componentNo;
  }
  
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << nValues << " entries.";
  values.resize(nValues);
  VecGetValues(this->values_, nValues, indices.data(), values.data());
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  // store the array indices for values_ array in dofGlobalNo
  for (int i=0; i<N; i++)
  {
    dofGlobalNo[i] = dofGlobalNo[i]*nComponents + componentNo;
  }
  
  VecGetValues(this->values_, N, (PetscInt *)dofGlobalNo.data(), values.data());
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values)
{
  // store the array indices for values_ array in dofGlobalNo
  const int nValues = dofGlobalNo.size();
  for (int i=0; i<nValues; i++)
  {
    dofGlobalNo[i] = dofGlobalNo[i]*nComponents + componentNo;
  }
  
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VecGetValues(this->values_, nValues, (PetscInt *)dofGlobalNo.data(), values.data()+previousSize);
}

template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  std::array<int,N*nComponents> indices;
  std::array<double,N*nComponents> result;
  
  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int i=0; i<N; i++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      indices[j] = dofGlobalNo[i]*nComponents + componentIndex;
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
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  std::array<PetscInt,nDofsPerElement> indices;
  
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    indices[dofIndex] = this->mesh_->getDofNo(elementNo,dofIndex)*nComponents + componentNo;
  }
  
  VecGetValues(this->values_, nDofsPerElement, indices.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  std::array<int,nDofsPerElement*nComponents> indices;
  std::array<double,nDofsPerElement*nComponents> result;
  
  //VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents;
  
  // prepare lookup indices for PETSc vector values_
  int j=0;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      indices[j] = this->mesh_->getDofNo(elementNo,dofIndex)*nComponents + componentIndex;
    }
  }
  
  VecGetValues(this->values_, nDofsPerElement*nComponents, indices.data(), result.data());
  
  //VLOG(2) << " indices: " << indices << ", retrieved values: " << result;
  
  // copy result to output values
  for (int dofIndex=0; dofIndex<nDofsPerElement; dofIndex++)
  {
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++, j++)
    {
      values[dofIndex][componentIndex] = result[dofIndex*nComponents+componentIndex];
      //VLOG(2) << "getElementValues element " << elementNo << ", dofIndex " << dofIndex << " componentIndex " << componentIndex << " value: " << values[dofIndex][componentIndex];
    }
  }
}

//! for a specific component, get a single value from global dof no.
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofGlobalNo)
{
  PetscInt index = dofGlobalNo*nComponents + componentNo;
  
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
  
  // loop over dof numbers
  int i=0;
  for (std::vector<dof_no_t>::iterator iter = dofGlobalNos.begin(); iter != dofGlobalNos.end(); iter++, i++)
  {  
    dof_no_t dofGlobalNo = *iter;
    
    // prepare lookup indices for PETSc vector values_
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = dofGlobalNo*nComponents + componentIndex;
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
  
  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    indices[componentIndex] = dofGlobalNo*nComponents + componentIndex;
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
