#include "field_variable/structured/01_field_variable_data_structured.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>
#include <array>

#include "mesh/structured.h"

namespace FieldVariable
{
 /*
template<typename BasisOnMeshType>
template<typename FieldVariableType>
void FieldVariableDataStructured<BasisOnMeshType>::
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
FieldVariableDataStructured<BasisOnMeshType>::
FieldVariableDataStructured() : FieldVariableBase<BasisOnMeshType>::FieldVariableBase()
{
}
  
//! contructor as data copy with a different name (component names are the same)
template<typename BasisOnMeshType>
FieldVariableDataStructured<BasisOnMeshType>::
FieldVariableDataStructured(FieldVariable<BasisOnMeshType> &rhs, std::string name) : 
  FieldVariableBase<BasisOnMeshType>::FieldVariableBase()
{
  // initialize everything from other field variable
  initializeFromFieldVariable(rhs, name, rhs.componentNames());
  
  // copy entries in values vector
  VecCopy(rhs.values(), this->values_);
}
  
//! constructor with mesh, name and components
template<typename BasisOnMeshType>
FieldVariableDataStructured<BasisOnMeshType>::
FieldVariableDataStructured(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames) : 
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
  
  
  LOG(DEBUG) << "FieldVariableDataStructured constructor, name=" << this->name_
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
FieldVariableDataStructured<BasisOnMeshType>::
~FieldVariableDataStructured()
{
  if(this->values_)
  {
    //PetscErrorCode ierr = VecDestroy(&this->values_); CHKERRV(ierr);
  }
}
  
template<typename BasisOnMeshType>
int FieldVariableDataStructured<BasisOnMeshType>::
nComponents() const
{
  return this->nComponents_;
}

template<typename BasisOnMeshType>
std::vector<std::string> FieldVariableDataStructured<BasisOnMeshType>::
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
std::array<element_no_t, BasisOnMeshType::Mesh::dim()> FieldVariableDataStructured<BasisOnMeshType>::
nElementsPerCoordinateDirection() const
{
  return this->mesh_->nElementsPerCoordinateDirection();
}

template<typename BasisOnMeshType>
std::size_t FieldVariableDataStructured<BasisOnMeshType>::
nEntries() const
{
  return this->nEntries_;
}

template<typename BasisOnMeshType>
dof_no_t FieldVariableDataStructured<BasisOnMeshType>::
nDofs() const
{
  assert(this->nComponents_ != 0);
  return this->nEntries_ / this->nComponents_;
}

template<typename BasisOnMeshType>
Vec &FieldVariableDataStructured<BasisOnMeshType>::
values()
{
  return this->values_;
}

template<typename BasisOnMeshType>
void FieldVariableDataStructured<BasisOnMeshType>::
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

//! write a exelem file header to a stream, for a particular element
template<typename BasisOnMeshType>
void FieldVariableDataStructured<BasisOnMeshType>::
outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  assert(0); // not implemented
}

//! write a exelem file header to a stream, for a particular element
template<typename BasisOnMeshType>
void FieldVariableDataStructured<BasisOnMeshType>::
outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
  assert(0); // not implemented
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<typename BasisOnMeshType>
bool FieldVariableDataStructured<BasisOnMeshType>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  return true;   // structured meshes do not have exfile representations 
}

template<typename BasisOnMeshType>
void FieldVariableDataStructured<BasisOnMeshType>::
output(std::ostream &stream) const
{
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::map<std::string, int> componentIndex_;   ///< names of the components and the component index (numbering starts with 0)
  int nComponents_;    ///< number of components
  std::size_t nEntries_;       ///< number of entries the PETSc vector values_ will have (if it is used). This number of dofs * nComponents
  
 
  stream << "\"" << this->name_ << "\", nEntries: " << nEntries_ 
    << ", isGeometryField: " << std::boolalpha << isGeometryField_ << std::endl
    << "  components:" << std::endl;
  for (auto &component : this->component_)
  {
    stream << component.second;
  }
}

};
