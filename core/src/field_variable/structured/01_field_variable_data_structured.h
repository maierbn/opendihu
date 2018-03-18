#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>
#include <petscvec.h>
#include "easylogging++.h"

#include "field_variable/00_field_variable_base.h"

namespace FieldVariable
{
 
// forward declaration
template<typename BasisOnMeshType> class FieldVariable;

/** Field variable for a structured mesh, i.e. dof and node information are purely implicit.
 *  This is used for RegularFixed and StructuredDeformable meshes.
 */
template<typename BasisOnMeshType>
class FieldVariableDataStructured : 
  public FieldVariableBase<BasisOnMeshType>
{
public:
  //! inherited constructor 
  //using FieldVariableBase<BasisOnMeshType>::FieldVariableBase;
 
  //! empty contructor
  FieldVariableDataStructured();
  
  //! contructor as data copy with a different name (component names are the same)
  FieldVariableDataStructured(FieldVariable<BasisOnMeshType> &rhs, std::string name);
  
  //! constructor with mesh, name and components
  FieldVariableDataStructured(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames);
 
  //! destructor
  virtual ~FieldVariableDataStructured();
 
  //! set all data but the values from a second field variable
  template<typename FieldVariableType>
  void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
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
  
  //! get the number of components
  int nComponents() const;
  
  //! get the names of the components
  std::vector<std::string> componentNames() const;
  
  //! get the number of elements per coordinate direction
  std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElementsPerCoordinateDirection() const;
  
  //! write a exelem file header to a stream, for a particular element, fieldVariableNo is the field index x) in the exelem file header
  void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo=-1);

  //! write a exelem file header to a stream, for a particular node
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2);

  //! get the internal PETSc vector values
  Vec &values();
  
  //! set all the data fields as well as the internal values PETSc vector
  void set(std::string name, std::vector<std::string> &componentNames, std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElements,
           std::size_t nEntries, bool isGeometryField, Vec &values);
  
  //! output string representation to stream for debugging
  void output(std::ostream &stream) const;
  
protected:
 
  //! get the number of entries of the internal values_ Vector
  std::size_t nEntries() const;
  
  //! get the number of dofs, i.e. the number of entries per component
  dof_no_t nDofs() const;
  
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::map<std::string, int> componentIndex_;   ///< names of the components and the component index (numbering starts with 0)
  int nComponents_;    ///< number of components
  std::size_t nEntries_;       ///< number of entries the PETSc vector values_ will have (if it is used). This number of dofs * nComponents
  
  Vec values_ = PETSC_NULL;          ///< Petsc vector containing the values, the values for the components are stored contiguous, e.g. (comp1val1, comp2val1, comp3val1, comp1val2, comp2val2, ...). Dof ordering proceeds fastest over dofs of a node, then over nodes, node numbering is along whole domain, fastes in x, then in y,z direction.
};
 
};  // namespace

#include "field_variable/structured/01_field_variable_data_structured.tpp"
