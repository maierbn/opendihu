#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/field_variable_base.h"
#include "field_variable/field_variable_interface.h"
#include "field_variable/field_variable_structured.h"
#include "field_variable/component.h"
#include "basis_on_mesh/04_basis_on_mesh_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/element_to_node_mapping.h"
#include "field_variable/node_to_dof_mapping.h"

namespace FieldVariable
{

/** FieldVariable class for UnstructuredDeformable mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> :
  public FieldVariableBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>,
  public Interface<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>
{ 
public:
  //! inherited constructor 
  using FieldVariableBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::FieldVariableBase;
  
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType> BasisOnMeshType;
  
  //! destructor
  virtual ~FieldVariable();
  
  //! set all data but the values from a second field variable
  template<typename FieldVariableType>
  void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames);
    
  //! set all data but the components from precomputed mappings
  void initializeFromMappings(std::string name, bool isGeometryField, 
                              std::shared_ptr<ExfileRepresentation> exfileRepresentation,
                              std::shared_ptr<ElementToDofMapping> elementToDofMapping,
                              std::shared_ptr<ElementToNodeMapping> elementToNodeMapping,
                              std::shared_ptr<NodeToDofMapping> nodeToDofMapping,
                              std::vector<std::string> componentNames);
  
  //! get the number of components
  int nComponents() const;
  
  //! get the number of elements
  element_idx_t nElements() const;
  
  //! get the number of nodes
  node_idx_t nNodes() const;
  
  //! get the names of the components
  std::vector<std::string> componentNames() const;
  
  //! get the exfile representation object
  std::shared_ptr<ExfileRepresentation> exfileRepresentation() const;
  
  //! get the element to dof mapping object
  std::shared_ptr<ElementToDofMapping> elementToDofMapping() const;
  
  //! get the element to dof mapping object
  std::shared_ptr<ElementToNodeMapping> elementToNodeMapping() const;  
  
  //! get the node to dof mapping object
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping() const;
  
  //! get the internal values vector
  Vec &values();
  
  //! get the number of scale factors
  int getNumberScaleFactors(element_idx_t globalElementNo) const;
  
  //! for a specific component, get values from their global dof no.s
  template<int N>
  void getValues(std::string component, std::array<int,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get values from their global dof no.s for all components
  template<int N, int nComponents>
  void getValues(std::array<int,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
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
  }
  
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(std::string component, element_idx_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
  {
    this->component_[component].getElementValues(elementNo, values);
  }
  
  //! get the values corresponding to all element-local dofs for all components
  template<int N, int nComponents>
  void getElementValues(element_idx_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
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
  
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_idx_t dofGlobalNo);

  //! get a single value from global dof no. for all components
  template<int nComponents>
  std::array<double,nComponents> getValue(node_idx_t dofGlobalNo);

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_idx_t currentElementGlobalNo, int fieldVariableNo=-1);
  
  //! write a exelem file header to a stream, for a particular node
  void outputHeaderExnode(std::ostream &file, node_idx_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1);
  
  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2);
  
  friend class BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>;
  
  //! resize internal representation variable to number of elements
  void setNumberElements(element_idx_t nElements);
  
  //! parse current component's exfile representation from file contents
  void parseHeaderFromExelemFile(std::string content);
  
  //! parse single element from exelem file
  void parseElementFromExelemFile(std::string content);
  
  //! read in values frorm exnode file
  void parseFromExnodeFile(std::string content);
  
  //! reduce memory consumption by removing duplicates in ExfileRepresentations
  void unifyMappings(std::shared_ptr<::FieldVariable::ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode);
  
  //! eliminate duplicate elementToDof and exfileRepresentation objects in components of two field variables (this and one other)
  void unifyMappings(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> &fieldVariable);
  
  //! initialize PETSc vector with size of total number of dofs for all components of this field variable
  void initializeValuesVector();
  
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  int getDofNo(element_idx_t elementNo, int dofIndex) const;
  
  //! return the component
  Component<BasisOnMeshType> &component(std::string key); 
  
  //! return the map of components
  std::map<std::string, Component<BasisOnMeshType>> &component();
  
  //! get the number of entries of the internal values_ Vector
  std::size_t nEntries() const;
  
  //! get the number of dofs, i.e. the number of entries per component
  dof_idx_t nDofs() const;
  
  //! multiply dof values with scale factors such that scale factor information is completely contained in dof values
  void eliminateScaleFactors();
  
  //! if the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate"
  bool isGeometryField();
  
  //! output string representation to stream for debugging
  void output(std::ostream &stream) const;
  
protected:
 
  //! create the element to dof mapping at each component
  void createElementToDofMapping(std::shared_ptr<::FieldVariable::ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode);
  
  //! initialize components
  void initializeComponents(std::vector<std::string> &componentNames, std::string exfileBasisRepresentation);
  
  int exfileNo_;    ///< number of the fieldvariable in exelem file (index starts at 1)
  std::size_t nEntries_;       ///< number of entries
  element_idx_t nElements_;    ///< number of elements
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::map<std::string, Component<BasisOnMeshType>> component_;    ///< one or multiple components of which this field variable consists of, with their name as key (e.g. 'x','y','z')
  std::shared_ptr<ExfileRepresentation> exfileRepresentation_;       ///< the indexing given in the exelem file, this is the same for all components
  std::shared_ptr<ElementToDofMapping> elementToDofMapping_;       ///< the element to dof mapping of all components, this is the same for all components
  std::shared_ptr<ElementToNodeMapping> elementToNodeMapping_;      ///< mapping from element-local node indices to global node numbers
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping_;       ///< the node to dof mapping of all components, this is the same for all components
  std::shared_ptr<Vec> values_;     ///< the vector that contains all values, the entries of all components are interleaved, e.g. (val1comp1, val1comp2, val2comp1, val2comp2, ...)
};

// output operator
template<int D, typename BasisFunctionType>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> &rhs);

};  // namespace

#include "field_variable/field_variable_unstructured_deformable.tpp"
