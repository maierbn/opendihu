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
  using FieldVariableBase<BasisOnMeshType>::FieldVariableBase;
 
  //! destructor
  virtual ~FieldVariableStructured();
 
  //! set all data but the values from a second field variable
  template<typename FieldVariableType>
  void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames);
  
  //! get the number of components
  int nComponents() const;
  
  //! get the names of the components
  std::vector<std::string> componentNames() const;
  
  //! get the number of elements
  std::array<int, BasisOnMeshType::Mesh::dim()> nElementsPerDimension() const;
  
  //! for a specific component, get values from their global dof no.s
  template<int N>
  void getValues(std::string component, std::array<int,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get values from their global dof no.s for all components
  template<int N, int nComponents>
  void getValues(std::array<int,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
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
  void getElementValues(std::string component, element_idx_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get the values corresponding to all element-local dofs for all components
  template<int nComponents>
  void getElementValues(element_idx_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
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
        indices[j] = BasisOnMeshType::getDofNo(this->nElements_,elementNo,dofIndex)*nComponents + componentIndex;
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
  double getValue(std::string component, node_idx_t dofGlobalNo);

  //! get a single value from global dof no. for all components
  template<int nComponents>
  std::array<double,nComponents> getValue(node_idx_t dofGlobalNo);

  //! write a exelem file header to a stream, for a particular element, fieldVariableNo is the field index x) in the exelem file header
  void outputHeaderExelem(std::ostream &file, element_idx_t currentElementGlobalNo, int fieldVariableNo=-1);

  //! write a exelem file header to a stream, for a particular node
  void outputHeaderExnode(std::ostream &file, node_idx_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2);

  //! get the internal PETSc vector values
  Vec &values();
  
  //! set all the data fields as well as the internal values PETSc vector
  void set(std::string name, std::vector<std::string> &componentNames, std::array<int, BasisOnMeshType::Mesh::dim()> nElements,
           int nEntries, bool isGeometryField, Vec &values);
  
protected:
 
  //! get the number of entries of the internal values_ Vector
  int nEntries() const;
  
  //! get the number of dofs, i.e. the number of entries per component
  int nDofs() const;
  
  std::array<int, BasisOnMeshType::Mesh::dim()> nElements_;    ///< number of elements in each coordinate direction
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::map<std::string, int> componentIndex_;   ///< names of the components and the component index (numbering starts with 0)
  int nComponents_;    ///< number of components
  int nEntries_;       ///< number of entries the PETSc vector values_ will have (if it is used). This number of dofs * nComponents
  
  Vec values_ = PETSC_NULL;          ///< Petsc vector containing the values, the values for the components are stored contiguous, e.g. (comp1val1, comp2val1, comp3val1, comp1val2, comp2val2, ...). Dof ordering proceeds fastest over dofs of a node, then over nodes, node numbering is along whole domain, fastes in x, then in y,z direction.
};
 
};  // namespace

#include "field_variable/field_variable_structured.tpp"
