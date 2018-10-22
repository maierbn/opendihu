#include "field_variable/structured/02_field_variable_data_structured.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>
#include <array>

#include "mesh/structured.h"

namespace FieldVariable
{
/*
template<typename FunctionSpaceType, int nComponents>
FieldVariableDataStructured<FunctionSpaceType,nComponents>::
FieldVariableDataStructured() :
  FieldVariableComponents<FunctionSpaceType,nComponents>::FieldVariableComponents()
{
}*/

//! contructor as data copy with a different name (component names are the same)
template<typename FunctionSpaceType, int nComponents>
FieldVariableDataStructured<FunctionSpaceType,nComponents>::
FieldVariableDataStructured(FieldVariable<FunctionSpaceType,nComponents> &rhs, std::string name) :
  FieldVariableComponents<FunctionSpaceType,nComponents>::FieldVariableComponents()
{
  // initialize everything from other field variable
  this->componentNames_.resize(rhs.componentNames().size());
  std::copy(rhs.componentNames().begin(), rhs.componentNames().end(), this->componentNames_.begin());
  
  this->name_ = name;
  this->isGeometryField_ = false;
  this->functionSpace_ = rhs.functionSpace();

  assert(this->functionSpace_);
  

  // create new distributed petsc vec as copy of rhs values vector
  if (rhs.partitionedPetscVec())
  {
    // if rhs is not a geometry field an therefore has a partitionedPetscVec, use that
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpaceType,nComponents>>(*rhs.partitionedPetscVec(), name);
  }
  else
  {
    // else create new values vector
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpaceType,nComponents>>(this->functionSpace_->meshPartition(), name);
  }
}

//! contructor as data copy with a different name and different components
template<typename FunctionSpaceType, int nComponents>
template <int nComponents2>
FieldVariableDataStructured<FunctionSpaceType,nComponents>::
FieldVariableDataStructured(FieldVariable<FunctionSpaceType,nComponents2> &rhs, std::string name, std::vector<std::string> componentNames) :
  FieldVariableComponents<FunctionSpaceType,nComponents>::FieldVariableComponents()
{
  // initialize everything from other field variable
  assert(componentNames.size() == nComponents);
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());
  
  this->name_ = name;
  this->isGeometryField_ = false;
  this->functionSpace_ = rhs.functionSpace();

  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());
  
  // create new distributed petsc vec as copy of rhs values vector
  if (rhs.partitionedPetscVec())
  {
    // if rhs is not a geometry field an therefore has a partitionedPetscVec, use that
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpaceType,nComponents>>(*rhs.partitionedPetscVec(), name);
  }
  else
  {
    // else create new values vector
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpaceType,nComponents>>(this->functionSpace_->meshPartition(), name);
  }
}

//! constructor with functionSpace, name and components
template<typename FunctionSpaceType, int nComponents>
FieldVariableDataStructured<FunctionSpaceType,nComponents>::
FieldVariableDataStructured(std::shared_ptr<FunctionSpaceType> functionSpace, std::string name, std::vector<std::string> componentNames, bool isGeometryField) :
  FieldVariableComponents<FunctionSpaceType,nComponents>::FieldVariableComponents()
{
  this->name_ = name;
  this->isGeometryField_ = isGeometryField;
  this->functionSpace_ = functionSpace;

  assert(this->functionSpace_);
  
  //std::shared_ptr<Mesh::Structured<FunctionSpaceType::dim()>> meshStructured = std::static_pointer_cast<Mesh::Structured<FunctionSpaceType::dim()>>(functionSpace);

  // assign component names
  assert(nComponents == componentNames.size());
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());

  LOG(DEBUG) << "FieldVariableDataStructured constructor, name=" << this->name_
   << ", components: " << nComponents << ", isGeometryField: " << this->isGeometryField_;

  bool isStructuredRegularFixed = std::is_same<typename FunctionSpaceType::Mesh, Mesh::StructuredRegularFixedOfDimension<1>>::value
    || std::is_same<typename FunctionSpaceType::Mesh, Mesh::StructuredRegularFixedOfDimension<2>>::value
    || std::is_same<typename FunctionSpaceType::Mesh, Mesh::StructuredRegularFixedOfDimension<3>>::value;

  // create a new values vector for the new field variable
  if (!this->isGeometryField_ || !isStructuredRegularFixed)
  {
    LOG(DEBUG) << "create a field variable with values_ vector";
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpaceType,nComponents>>(this->functionSpace_->meshPartition(), name);
  }
  else
  {
    LOG(DEBUG) << "create a geometry field variable without values_ vector (because it isStructuredRegularFixed)";
  }
}

template<typename FunctionSpaceType, int nComponents>
FieldVariableDataStructured<FunctionSpaceType,nComponents>::
~FieldVariableDataStructured()
{
}
/*
template<typename FunctionSpaceType, int nComponents>
std::array<element_no_t, FunctionSpaceType::Mesh::dim()> FieldVariableDataStructured<FunctionSpaceType,nComponents>::
nElementsPerCoordinateDirectionLocal() const
{
  return this->functionSpace_->nElementsPerCoordinateDirectionLocal();
}*/
/*
template<typename FunctionSpaceType, int nComponents>
dof_no_t FieldVariableDataStructured<FunctionSpaceType,nComponents>::
nDofsLocalWithoutGhosts() const
{
  // this is the same as this->nEntries_ / nComponents only if it is not the geometry mesh of StructuredRegularFixed mesh
  return this->functionSpace_->nDofsLocalWithoutGhosts();
}*/

template<typename FunctionSpaceType, int nComponents>
Vec &FieldVariableDataStructured<FunctionSpaceType,nComponents>::
valuesLocal(int componentNo)
{
  assert(this->values_);
  return this->values_->valuesLocal(componentNo);
}

template<typename FunctionSpaceType, int nComponents>
Vec &FieldVariableDataStructured<FunctionSpaceType,nComponents>::
valuesGlobal(int componentNo)
{
  assert(this->values_);
  return this->values_->valuesGlobal(componentNo);
}

template<typename FunctionSpaceType, int nComponents>
Vec &FieldVariableDataStructured<FunctionSpaceType,nComponents>::
getValuesContiguous()
{
  assert(this->values_);
  return this->values_->getValuesContiguous();
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableDataStructured<FunctionSpaceType,nComponents>::
restoreValuesContiguous()
{
  assert(this->values_);
  this->values_->restoreValuesContiguous();
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents>> FieldVariableDataStructured<FunctionSpaceType,nComponents>::
partitionedPetscVec()
{
  return this->values_; 
}


/*
 * why is this needed?
template<typename FunctionSpaceType, int nComponents>
void FieldVariableDataStructured<FunctionSpaceType,nComponents>::
initialize(std::string name, std::vector<std::string> &componentNames, bool isGeometryField)
{
  // set member variables from arguments
  this->name_ = name;
  this->isGeometryField_ = isGeometryField;

  // copy component names
  assert(nComponents == (int)componentNames.size());
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());

  // if there is not yet a values vector, create PartitionedPetscVec (not for RegularFixed mesh geometry fields)
  if (this->values_ == nullptr)
  {
    initializeValuesVector();
  }
}
*/

//! write a exelem file header to a stream, for a particular element
template<typename FunctionSpaceType, int nComponents>
void FieldVariableDataStructured<FunctionSpaceType,nComponents>::
outputHeaderExelem(std::ostream &stream, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  // output first line of header
  stream << " " << (fieldVariableNo+1) << ") " << this->name_ << ", " << (this->isGeometryField_? "coordinate" : "field")
    << ", rectangular cartesian, #Components=" << nComponents << std::endl;

  // output headers of components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::string name = this->componentNames_[componentNo];

    // compose basis representation string, e.g. "l.Lagrange*l.Lagrange"
    std::stringstream basisFunction;
    switch(FunctionSpaceType::BasisFunction::getBasisOrder())
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
    basisFunction << FunctionSpaceType::BasisFunction::getBasisFunctionString();

    stream << " " << name << ".   " << StringUtility::multiply<FunctionSpaceType::dim()>(basisFunction.str())
      << ", no modify, standard node based." << std::endl
      << "   #Nodes= " << FunctionSpaceType::nNodesPerElement() << std::endl;

    // loop over nodes of a representative element
    for (int nodeIndex = 0; nodeIndex < FunctionSpaceType::nNodesPerElement(); nodeIndex++)
    {
      stream << "   " << (nodeIndex+1) << ".  #Values=" << FunctionSpaceType::nDofsPerNode() << std::endl
        << "      Value indices: ";
      
      for (int dofIndex = 0; dofIndex <FunctionSpaceType::nDofsPerNode(); dofIndex++)
      {
        stream << " " << dofIndex + 1;
      }
      
      stream << "      Scale factor indices: 0"  << std::endl;
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
template<typename FunctionSpaceType, int nComponents>
void FieldVariableDataStructured<FunctionSpaceType,nComponents>::
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
  stream << " " << (fieldVariableNo+1) << ") " << this->name_ << ", " << (this->isGeometryField_? "coordinate" : "field")
    << ", rectangular cartesian, #Components=" << nComponents << std::endl;

  // output headers of components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::string name = this->componentNames_[componentNo];
    stream << "  " << name << ".  Value index=" << valueIndex+1 << ", #Derivatives= 0" << std::endl;

    valueIndex += FunctionSpaceType::nDofsPerNode();
  }
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<typename FunctionSpaceType, int nComponents>
bool FieldVariableDataStructured<FunctionSpaceType,nComponents>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  return true;   // structured meshes do not have exfile representations
}

template<typename FunctionSpaceType, int nComponents>
void FieldVariableDataStructured<FunctionSpaceType,nComponents>::
output(std::ostream &stream) const
{
  // only output if on rank 0
  if (this->functionSpace_->meshPartition()->ownRankNo() == 0)
  {
    stream << "\"" << this->name_ << "\""
      << ", isGeometryField: " << std::boolalpha << this->isGeometryField_
      << ", " << this->componentNames_.size() << (this->componentNames_.size() == 1? " component:" : " components: ");
    for (auto &componentName : this->componentNames_)
    {
      stream << "\"" << componentName << "\", ";
    }
  }

  if (values_)
  {
    // the values have to be output from all ranks
    stream << *values_;
  }
  else
  {
    if (this->functionSpace_->meshPartition()->ownRankNo() == 0)
    {
      stream << "(values not set)";
    }
  }
}

};
