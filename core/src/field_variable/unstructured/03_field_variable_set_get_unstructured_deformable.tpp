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

//! get values from their global dof no.s for all components, this eventually does not get all values if there are multiple versions
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  // get number of dofs
  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nLocalDofs();

  std::array<double,nComponents> resultVector;

  // transform global dof no.s to vector indices of first component
  for (int valueIndex = 0; valueIndex < N; valueIndex++)
  {
    // TODO: map to global no
    int valuesVectorIndex = dofGlobalNo[valueIndex]*nComponents;

    // create indices vector
    std::array<int,nComponents> indices;
    for(int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = componentIndex*nDofs + dofGlobalNo[valueIndex];
    }
    
    // get values and assign them to result values vector
    VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
    values[valueIndex] = resultVector;
  }
}

//! for a specific component, get all values
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getValues(values, onlyNodalValues);
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].template getValues<N>(dofGlobalNo, values);
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getValues(dofGlobalNo, values);
}

//! for a specific component, get a single value from global dof no.
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofGlobalNo)
{
  assert(componentNo >= 0 && componentNo < nComponents);
  return this->component_[componentNo].getValue(dofGlobalNo);
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nLocalElements());
  assert(componentNo >= 0 && componentNo < nComponents);
  
  this->component_[componentNo].getElementValues(elementNo, values);
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nLocalElements());
  
  // get number of dofs
  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nLocalDofs();

  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();

  // TODO: local to global
  const std::vector<dof_no_t> &dofGlobalNo = this->elementToDofMapping_->getElementDofs(elementNo);
  std::array<double,nComponents> resultVector;

  // transform global dof no.s to vector indices of first component
  for (int valueIndex = 0; valueIndex < nDofsPerElement; valueIndex++)
  {
    dof_no_t valuesVectorIndex = dofGlobalNo[valueIndex]*nComponents;

    // create indices vector
    std::array<PetscInt,nComponents> indices;
    for(int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = componentIndex*nDofs + dofGlobalNo[valueIndex];
    }
    
    // get values and assign them to result values vector
    VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
    values[valueIndex] = resultVector;
  }
}

//! copy the values from another field variable of the same type
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  VecCopy(*rhs.values_, *this->values_);
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->nComponents == 1);
  const int nValues = values.size();

  VecSetValues(*this->values_, nValues, (const int *) dofGlobalNos.data(), values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

//! set value for all dofs
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(double value)
{
  VecSet(*this->values_, value);
}

//! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode)
{
  // get number of dofs
  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nLocalDofs();

  std::array<int,nComponents> indices;

  // loop over dof numbers
  int i=0;
  for (std::vector<dof_no_t>::iterator iter = dofGlobalNos.begin(); iter != dofGlobalNos.end(); iter++, i++)
  {
    dof_no_t dofGlobalNo = *iter;

    // prepare lookup indices for PETSc vector values_
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      indices[componentNo] = componentNo*nDofs + dofGlobalNo;
    }

    VecSetValues(*this->values_, nComponents, indices.data(), values[i].data(), petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

//! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode)
{
  // get number of dofs
  assert(this->mesh_);
  const dof_no_t nDofs = this->mesh_->nLocalDofs();

  std::array<int,nComponents> indices;

  // prepare lookup indices for PETSc vector values_
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    indices[componentNo] = componentNo*nDofs + dofGlobalNo;
  }

  VecSetValues(*this->values_, nComponents, indices.data(), value.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

//! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetUnstructured<BasisOnMeshType,nComponents>::
flushSetValues()
{
  VecAssemblyBegin(*this->values_);
  VecAssemblyEnd(*this->values_);
}

};  // namespace
