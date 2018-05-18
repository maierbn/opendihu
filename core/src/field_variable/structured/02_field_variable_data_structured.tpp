#include "field_variable/structured/02_field_variable_data_structured.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>
#include <array>

#include "mesh/structured.h"

namespace FieldVariable
{

template<typename BasisOnMeshType, int nComponents>
FieldVariableDataStructured<BasisOnMeshType,nComponents>::
FieldVariableDataStructured() :
  FieldVariableComponents<BasisOnMeshType,nComponents>::FieldVariableComponents()
{
}

//! contructor as data copy with a different name (component names are the same)
template<typename BasisOnMeshType, int nComponents>
FieldVariableDataStructured<BasisOnMeshType,nComponents>::
FieldVariableDataStructured(FieldVariable<BasisOnMeshType,nComponents> &rhs, std::string name) :
  FieldVariableComponents<BasisOnMeshType,nComponents>::FieldVariableComponents()
{
  // initialize everything from other field variable
  std::vector<std::string> componentNames(rhs.componentNames().size());
  std::copy(rhs.componentNames().begin(), rhs.componentNames().end(), componentNames.begin());

  initializeFromFieldVariable(rhs, name, componentNames);

  // copy entries in values vector
  VecCopy(rhs.values(), this->values_);
}

//! constructor with mesh, name and components
template<typename BasisOnMeshType, int nComponents>
FieldVariableDataStructured<BasisOnMeshType,nComponents>::
FieldVariableDataStructured(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames) :
  FieldVariableComponents<BasisOnMeshType,nComponents>::FieldVariableComponents()
{
  this->name_ = name;
  this->isGeometryField_ = false;
  this->mesh_ = mesh;

  std::shared_ptr<Mesh::Structured<BasisOnMeshType::dim()>> meshStructured = std::static_pointer_cast<Mesh::Structured<BasisOnMeshType::dim()>>(mesh);

  // assign component names
  assert(nComponents == componentNames.size());
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());

  this->nEntries_ = mesh->nDofs() * nComponents;

  LOG(DEBUG) << "FieldVariableDataStructured constructor, name=" << this->name_
   << ", components: " << nComponents << ", nEntries: " << this->nEntries_;

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

template<typename BasisOnMeshType, int nComponents>
FieldVariableDataStructured<BasisOnMeshType,nComponents>::
~FieldVariableDataStructured()
{
  if(this->values_)
  {
    //PetscErrorCode ierr = VecDestroy(&this->values_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType, int nComponents>
template<typename FieldVariableType>
void FieldVariableDataStructured<BasisOnMeshType,nComponents>::
initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
{
  this->name_ = name;
  this->isGeometryField_ = false;
  this->mesh_ = fieldVariable.mesh();

  // copy component names
  assert(nComponents == (int)componentNames.size());
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());

  this->nEntries_ = fieldVariable.nDofs() * nComponents;

  LOG(DEBUG) << "FieldVariable::initializeFromFieldVariable, name=" << this->name_
   << ", components: " << nComponents << ", Vec nEntries: " << this->nEntries_;

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

template<typename BasisOnMeshType, int nComponents>
std::array<element_no_t, BasisOnMeshType::Mesh::dim()> FieldVariableDataStructured<BasisOnMeshType,nComponents>::
nElementsPerCoordinateDirection() const
{
  return this->mesh_->nElementsPerCoordinateDirection();
}

template<typename BasisOnMeshType, int nComponents>
std::size_t FieldVariableDataStructured<BasisOnMeshType,nComponents>::
nEntries() const
{
  return this->nEntries_;
}

template<typename BasisOnMeshType, int nComponents>
dof_no_t FieldVariableDataStructured<BasisOnMeshType,nComponents>::
nDofs() const
{
  return this->nEntries_ / nComponents;
}

template<typename BasisOnMeshType, int nComponents>
Vec &FieldVariableDataStructured<BasisOnMeshType,nComponents>::
values()
{
  return this->values_;
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableDataStructured<BasisOnMeshType,nComponents>::
set(std::string name, std::vector<std::string> &componentNames, std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElements,
    std::size_t nEntries, bool isGeometryField, Vec &values)
{
  this->name_ = name;
  this->isGeometryField_ = isGeometryField;

  // copy component names
  assert(nComponents == (int)componentNames.size());
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());

  nEntries_ = nEntries;
  values_ = values;
}

//! write a exelem file header to a stream, for a particular element
template<typename BasisOnMeshType, int nComponents>
void FieldVariableDataStructured<BasisOnMeshType,nComponents>::
outputHeaderExelem(std::ostream &stream, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  // output first line of header
  stream << " " << (fieldVariableNo+1) << ") " << this->name_ << ", " << (isGeometryField_? "coordinate" : "field")
    << ", rectangular cartesian, #Components=" << nComponents << std::endl;

  // output headers of components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::string name = this->componentNames_[componentNo];

    // compose basis representation string, e.g. "l.Lagrange*l.Lagrange"
    std::stringstream basisFunction;
    switch(BasisOnMeshType::BasisFunction::getBasisOrder())
    {
    case 1:
      basisFunction << "l.";
      break;
    case 2:
      basisFunction << "q.";
      break;
    case 3:
      basisFunction << "c.";
      break;
    default:
      break;
    }
    basisFunction << BasisOnMeshType::BasisFunction::getBasisFunctionString();

    stream << " " << name << ".   " << StringUtility::multiply<BasisOnMeshType::dim()>(basisFunction.str())
      << ", no modify, standard node based." << std::endl;

    // loop over nodes of a representative element
    for (int nodeIndex = 0; nodeIndex < BasisOnMeshType::nNodesPerElement(); nodeIndex++)
    {
      stream << "   " << (nodeIndex+1) << ".  #Values=" << BasisOnMeshType::nDofsPerNode() << std::endl
        << "      Scale factor indices: 0"  << std::endl;
    }
  }
  /*
   1) Geometry, coordinate, rectangular cartesian, #Components=2
     x.   l.Lagrange*l.Lagrange, no modify, standard node based.
     #Nodes= 4
     1.  #Values=1
      Value indices:     1
      Scale factor indices:    1
     2.  #Values=1
      Value indices:     1
      Scale factor indices:    2
   */
}

//! write a exnode file header to a stream, for a particular element
template<typename BasisOnMeshType, int nComponents>
void FieldVariableDataStructured<BasisOnMeshType,nComponents>::
outputHeaderExnode(std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
/* example output:
 1) Geometry, coordinate, rectangular cartesian, #Components=2
   x.  Value index= 1, #Derivatives= 0
   y.  Value index= 2, #Derivatives= 0
 2) ScalarField, field,  rectangular cartesian, #Components=1
   1.  Value index= 3, #Derivatives= 0
*/
  // output first line of header
  stream << " " << (fieldVariableNo+1) << ") " << this->name_ << ", " << (isGeometryField_? "coordinate" : "field")
    << ", rectangular cartesian, #Components=" << nComponents << std::endl;

  // output headers of components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::string name = this->componentNames_[componentNo];
    stream << "  " << name << ".  Value index=" << valueIndex+1 << ", #Derivatives= 0" << std::endl;

    valueIndex += BasisOnMeshType::nDofsPerNode();
  }
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<typename BasisOnMeshType, int nComponents>
bool FieldVariableDataStructured<BasisOnMeshType,nComponents>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  return true;   // structured meshes do not have exfile representations
}

template<typename BasisOnMeshType, int nComponents>
void FieldVariableDataStructured<BasisOnMeshType,nComponents>::
output(std::ostream &stream) const
{
  stream << "\"" << this->name_ << "\", nEntries: " << nEntries_
    << ", isGeometryField: " << std::boolalpha << isGeometryField_ << std::endl
    << "  components:" << std::endl;
  for (auto &componentName : this->componentNames_)
  {
    stream << componentName;
  }
}

//! if the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate"
template<typename BasisOnMeshType, int nComponents>
bool FieldVariableDataStructured<BasisOnMeshType,nComponents>::isGeometryField() const
{
  return isGeometryField_;
}

};
