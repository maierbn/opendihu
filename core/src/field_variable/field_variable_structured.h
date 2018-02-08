#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/field_variable_base.h"
#include "field_variable/field_variable_interface.h"
#include "field_variable/component.h"
#include "basis_on_mesh/04_basis_on_mesh_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/element_to_node_mapping.h"
#include "field_variable/node_to_dof_mapping.h"

namespace FieldVariable
{

/** Field variable for a structured mesh, i.e. dof and node information are purely implicit.
 *  This is used for RegularFixed and StructuredDeformable meshes.
 */
template<typename BasisOnMeshType>
class FieldVariableStructured : public FieldVariableBase<BasisOnMeshType>
{
public:
  //! inherited constructor 
  //using FieldVariableBase<BasisOnMeshType>::FieldVariableBase;
 
  //! empty contructor
  FieldVariableStructured();
  
  //! contructor as data copy with a different name (component names are the same)
  FieldVariableStructured(FieldVariable<BasisOnMeshType> &rhs, std::string name);
  
  //! constructor with mesh, name and components
  FieldVariableStructured(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames);
 
  //! destructor
  virtual ~FieldVariableStructured();
 
  //! set all data but the values from a second field variable
  template<typename FieldVariableType>
  void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
  {
    this->name_ = name;
    this->nElementsPerDimension_ = fieldVariable.nElementsPerDimension();
    this->isGeometryField_ = false;
    this->mesh_ = fieldVariable.mesh();
    
    int index = 0;
    for (auto &componentName : componentNames)
    {
      this->componentIndex_.insert(std::pair<std::string,int>(componentName,index++));
    }
    this->nComponents_ = componentNames.size();
    this->nEntries_ = fieldVariable.nDofs() * this->nComponents_;
    
    
    LOG(DEBUG) << "FieldVariable::initializeFromFieldVariable, name=" << this->name_ << ", nElements: " << this->nElementsPerDimension_
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
  
  //! get the number of elements
  std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElementsPerDimension() const;
  
  //! for a specific component, get all values
  void getValues(std::string component, std::vector<double> &values);
  
  //! for a specific component, get values from their global dof no.s
  template<int N>
  void getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get values from their global dof no.s for all components
  template<int N, int nComponents>
  void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
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
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(std::string component, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get the values corresponding to all element-local dofs for all components
  template<int nComponents>
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
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
        indices[j] = BasisOnMeshType::getDofNo(this->nElementsPerDimension_,elementNo,dofIndex)*nComponents + componentIndex;
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
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_no_t dofGlobalNo);

  //! get a single value from global dof no. for all components
  template<int nComponents>
  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);

  //! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<int nComponents>
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values)
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
  }

  //! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<int nComponents>
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value)
  {
    std::array<int,nComponents> indices;
    
    // prepare lookup indices for PETSc vector values_
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = dofGlobalNo*this->nComponents_ + componentIndex;
    }
    
    VecSetValues(this->values_, nComponents, indices.data(), value.data(), INSERT_VALUES);
    // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
  }

    
  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  void flushSetValues();
  
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
  
protected:
 
  //! get the number of entries of the internal values_ Vector
  std::size_t nEntries() const;
  
  //! get the number of dofs, i.e. the number of entries per component
  dof_no_t nDofs() const;
  
  std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElementsPerDimension_;    ///< number of elements in each coordinate direction
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::map<std::string, int> componentIndex_;   ///< names of the components and the component index (numbering starts with 0)
  int nComponents_;    ///< number of components
  std::size_t nEntries_;       ///< number of entries the PETSc vector values_ will have (if it is used). This number of dofs * nComponents
  
  Vec values_ = PETSC_NULL;          ///< Petsc vector containing the values, the values for the components are stored contiguous, e.g. (comp1val1, comp2val1, comp3val1, comp1val2, comp2val2, ...). Dof ordering proceeds fastest over dofs of a node, then over nodes, node numbering is along whole domain, fastes in x, then in y,z direction.
};
 
};  // namespace

#include "field_variable/field_variable_structured.tpp"
