#include "field_variable/field_variable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>
#include <array>

namespace FieldVariable
{
 /*
template<typename BasisOnMeshType>
template<typename FieldVariableType>
void FieldVariableStructured<BasisOnMeshType>::
initializeFromFieldVariable(const FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
{
  this->name_ = name;
  this->isGeometryField_ = false;
  this->mesh_ = fieldVariable.mesh();
  
  int index = 0;
  for (auto &componentName : componentNames)
  {
    this->componentIndex_.insert(std::pair<std::string,int>(componentName,index++));
  }
  this->nComponents_ = componentNames.size();
  this->nEntries_ = fieldVariable.nDofs() * this->nComponents_;
  
  
  LOG(DEBUG) << "FieldVariable::initializeFromFieldVariable, name=" << this->name_ 
   << ", components: " << this->nComponents_ << ", nEntries: " << this->nEntries_;
  
  assert(this->nEntries_ != 0);
   
  // create a new values vector for the new field variable
  
  // create vector
  PetscErrorCode ierr;
  // initialize PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &this->values_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) this->values_, this->name_.c_str()); CHKERRV(ierr);
  
  // initialize size of vector
  ierr = VecSetSizes(this->values_, PETSC_DECIDE, this->nEntries_); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(this->values_);  CHKERRV(ierr);
}
*/
template<typename BasisOnMeshType>
FieldVariableStructured<BasisOnMeshType>::
FieldVariableStructured() : FieldVariableBase<BasisOnMeshType>::FieldVariableBase()
{
}
  
//! contructor as data copy with a different name (component names are the same)
template<typename BasisOnMeshType>
FieldVariableStructured<BasisOnMeshType>::
FieldVariableStructured(FieldVariable<BasisOnMeshType> &rhs, std::string name) : 
  FieldVariableBase<BasisOnMeshType>::FieldVariableBase()
{
  // initialize everything from other field variable
  initializeFromFieldVariable(rhs, name, rhs.componentNames());
  
  // copy entries in values vector
  VecCopy(rhs.values(), this->values_);
}
  
//! constructor with mesh, name and components
template<typename BasisOnMeshType>
FieldVariableStructured<BasisOnMeshType>::
FieldVariableStructured(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames) : 
  FieldVariableBase<BasisOnMeshType>::FieldVariableBase()
{
  this->name_ = name;
  this->isGeometryField_ = false;
  this->mesh_ = mesh;
  
  std::shared_ptr<Mesh::Structured<BasisOnMeshType::dim()>> meshStructured = std::static_pointer_cast<Mesh::Structured<BasisOnMeshType::dim()>>(mesh);
  
  int index = 0;
  for (auto &componentName : componentNames)
  {
    this->componentIndex_.insert(std::pair<std::string,int>(componentName,index++));
  }
  this->nComponents_ = componentNames.size();
  this->nEntries_ = mesh->nDofs() * this->nComponents_;
  
  
  LOG(DEBUG) << "FieldVariableStructured constructor, name=" << this->name_
   << ", components: " << this->nComponents_ << ", nEntries: " << this->nEntries_;
  
  assert(this->nEntries_ != 0);
   
  // create a new values vector for the new field variable
  
  // create vector
  PetscErrorCode ierr;
  // initialize PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &this->values_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) this->values_, this->name_.c_str()); CHKERRV(ierr);
  
  // initialize size of vector
  ierr = VecSetSizes(this->values_, PETSC_DECIDE, this->nEntries_); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(this->values_);  CHKERRV(ierr);
}
 
template<typename BasisOnMeshType>
FieldVariableStructured<BasisOnMeshType>::
~FieldVariableStructured()
{
  if(this->values_)
  {
    //PetscErrorCode ierr = VecDestroy(&this->values_); CHKERRV(ierr);
  }
}
  
template<typename BasisOnMeshType>
int FieldVariableStructured<BasisOnMeshType>::
nComponents() const
{
  return this->nComponents_;
}

template<typename BasisOnMeshType>
std::vector<std::string> FieldVariableStructured<BasisOnMeshType>::
componentNames() const
{
  std::vector<std::string> result;
  result.reserve(this->nComponents_);
  for (auto componentIndex : this->componentIndex_)
  {   
    result.push_back(componentIndex.first);
  }
  return result;
}
  
template<typename BasisOnMeshType>
std::array<element_no_t, BasisOnMeshType::Mesh::dim()> FieldVariableStructured<BasisOnMeshType>::
nElementsPerCoordinateDirection() const
{
  return this->mesh_->nElementsPerCoordinateDirection();
}

template<typename BasisOnMeshType>
std::size_t FieldVariableStructured<BasisOnMeshType>::
nEntries() const
{
  return this->nEntries_;
}

template<typename BasisOnMeshType>
dof_no_t FieldVariableStructured<BasisOnMeshType>::
nDofs() const
{
  assert(this->nComponents_ != 0);
  return this->nEntries_ / this->nComponents_;
}

template<typename BasisOnMeshType>
Vec &FieldVariableStructured<BasisOnMeshType>::
values()
{
  return this->values_;
}

template<typename BasisOnMeshType>
void FieldVariableStructured<BasisOnMeshType>::
set(std::string name, std::vector<std::string> &componentNames, std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElements,
    std::size_t nEntries, bool isGeometryField, Vec &values)
{
  this->name_ = name;
  this->isGeometryField_ = isGeometryField;
  
  // create numbering for components
  int i = 0;
  for (auto& componentName : componentNames)
  {
    componentIndex_[componentName] = i++;
  }
  nComponents_ = componentNames.size();
  nEntries_ = nEntries;
  values_ = values;
}

//! for a specific component, get all values
template<typename BasisOnMeshType>
void FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
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
double FieldVariableStructured<BasisOnMeshType>::
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
std::array<double,nComponents> FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
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
void FieldVariableStructured<BasisOnMeshType>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->nComponents == 1);
  const int nValues = values.size();

  VecSetValues(this->values_, nValues, (const int *) dofGlobalNos.data(), values.data(), petscInsertMode);
  
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}

template<typename BasisOnMeshType>
void FieldVariableStructured<BasisOnMeshType>::
flushSetValues()
{
  VecAssemblyBegin(this->values_); 
  VecAssemblyEnd(this->values_);
}
 
//! write a exelem file header to a stream, for a particular element
template<typename BasisOnMeshType>
void FieldVariableStructured<BasisOnMeshType>::
outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  assert(0); // not implemented
}

//! write a exelem file header to a stream, for a particular element
template<typename BasisOnMeshType>
void FieldVariableStructured<BasisOnMeshType>::
outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
  assert(0); // not implemented
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<typename BasisOnMeshType>
bool FieldVariableStructured<BasisOnMeshType>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  return true;   // structured meshes do not have exfile representations 
}

};
