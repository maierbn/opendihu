#include "field_variable/unstructured/02_field_variable_set_get_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>
#include <iomanip>
#include <cassert>

namespace FieldVariable
{
  
using namespace StringUtility;

//! for a specific component, get all values
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValues(std::string component, std::vector<double> &values, bool onlyNodalValues)
{
  this->component_[component].getValues(values, onlyNodalValues);
}

//! for a specific component, get values from their global dof no.s
template<int D, typename BasisFunctionType>
template<int N>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  this->component_[component].getValues(dofGlobalNo, values);
}
/*
template<int D, typename BasisFunctionType>
template<int N, int nComponents>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  for (int valueIndex = 0; valueIndex < N; valueIndex++)
  {
    int valuesVectorIndex = dofGlobalNo[valueIndex]*nComponents;
    
    // create indices vector with values {0,1,2,...,nComponents-1}
    std::array<int,nComponents> indices;
    for(int i=0; i<nComponents; i++)
      indices[i] = valuesVectorIndex + i;
    
    // get values and assign them to result values vector
    VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
    values[valueIndex] = resultVector;
  }
}*/

/*
//! for a specific component, get the values corresponding to all element-local dofs
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getElementValues(std::string component, element_no_t elementNo, std::array<double,BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  this->component_[component].getElementValues(elementNo, values);
}
*/
/*
//! get the values corresponding to all element-local dofs for all components
template<int D, typename BasisFunctionType>
template<std::size_t nComponents>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  std::vector<int> &dofGlobalNo = this->elementToDofMapping_->getElementDofs(elementNo);
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  for (int valueIndex = 0; valueIndex < nDofsPerElement; valueIndex++)
  {
    int valuesVectorIndex = dofGlobalNo[valueIndex]*nComponents;
    
    // create indices vector with values {0,1,2,...,nComponents-1}
    std::array<int,nComponents> indices;
    for(int i=0; i<nComponents; i++)
      indices[i] = valuesVectorIndex + i;
    
    // get values and assign them to result values vector
    VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
    values[valueIndex] = resultVector;
  }
}
*/
//! for a specific component, get a single value from global dof no.
template<int D, typename BasisFunctionType>
double FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValue(std::string component, node_no_t dofGlobalNo)
{
  assert(this->component_.find(component) != this->component_.end());
  return this->component_[component].getValue(dofGlobalNo); 
}

//! get a single value from global dof no. for all components
template<int D, typename BasisFunctionType>
template<std::size_t nComponents>
std::array<double,nComponents> FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValue(node_no_t dofGlobalNo)
{
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  int valuesVectorIndex = dofGlobalNo*nComponents;
  
  // create indices vector with values {0,1,2,...,nComponents-1}
  std::array<int,nComponents> indices;
  for(int i=0; i<nComponents; i++)
    indices[i] = valuesVectorIndex + i;
  
  // get values and assign them to result values vector
  VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
  return resultVector;
}

//! copy the values from another field variable of the same type
template<int D,typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> &rhs)
{
  VecCopy(*rhs.values_, *this->values_);
}

/*
//! set values for dofs
template<int D,typename BasisFunctionType>
template<std::size_t nComponents>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
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
template<int D,typename BasisFunctionType>
template<std::size_t nComponents>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
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

template<int D,typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->nComponents == 1);
  const int nValues = values.size();

  VecSetValues(*this->values_, nValues, (const int *) dofGlobalNos.data(), values.data(), petscInsertMode);
  
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}

//! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
template<int D,typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::
flushSetValues()
{
  VecAssemblyBegin(*this->values_); 
  VecAssemblyEnd(*this->values_);
}

};  // namespace
