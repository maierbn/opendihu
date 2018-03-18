#pragma once

#include <Python.h>  // has to be the first included header

#include <petscmat.h>
#include <iostream>
#include <memory>
#include <array>

#include "field_variable/unstructured/element_to_dof_mapping.h"
#include "field_variable/unstructured/exfile_representation.h"

namespace FieldVariable
{

template<typename BasisOnMeshType>
class Component
{
public:
 
  //! initialize values
  void initialize(std::shared_ptr<Vec> values, int nComponents, int componentIndex, int nElements);
 
  //! set the internal values PETSc vector
  void setValuesVector(std::shared_ptr<Vec> values);
  
  //! parse current component's exfile representation from file contents
  void parseHeaderFromExelemFile(std::string content);
  
  //! parse a part of the exelem file that describes a single element
  void parseElementFromExelemFile(std::string content);
  
  //! assign values of all dofs for a node, valuesBegin is an iterator that can iterate over the needed number of values, the order is given by dofs numbering
  void setNodeValues(node_no_t nodeGlobalNo, std::vector<double>::iterator valuesBegin);
  
  //! assign values of all dofs for a node, the values are given by valuesBegin iterator as they appear in the exnode block, valuesBegin points to the beginning of the subblock for the particular component
  void setNodeValuesFromBlock(node_no_t nodeGlobalNo, std::vector<double>::iterator valuesBegin);
  
  //! assign the name of the component
  void setName(std::string name, std::string exfileBasisFunctionSpecification="");
  
  //! set the exfileRepresentation pointer
  void setExfileRepresentation(std::shared_ptr<ExfileRepresentation> exfileRepresentation);
  
  //! set the setelementToDofMapping and nodeToDofMapping pointers
  void setDofMappings(std::shared_ptr<ElementToDofMapping> elementToDofMapping, std::shared_ptr<NodeToDofMapping> nodeToDofMapping);
  
  //! get the exfileRepresentation object
  std::shared_ptr<ExfileRepresentation> exfileRepresentation();
  
  //! get the elementToDofMapping object
  std::shared_ptr<ElementToDofMapping> elementToDofMapping();
  
  //! get the nodeToDofMapping object
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping();
  
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the total number of dofs for this component
  dof_no_t nDofs() const;
  
  //! return the number of elements
  element_no_t nElements() const;
  
  //! return the exfileBasisFunctionSpecification
  std::string exfileBasisFunctionSpecification() const;
  
  //! get all values
  void getValues(std::vector<double> &values);
  
  //! get values from their global dof no.s
  template<int N>
  void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get the values corresponding to all element-local dofs
  void getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get a single value from global dof no.
  double getValue(node_no_t dofGlobalNo);
  
  //! get the number of scale factors for a given node
  int getNumberScaleFactors(node_no_t nodeGlobalNo) const;
  
  //! write an exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo);
  
  //! write an exnode file header to a stream, for a particular node
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex);
  
  //! output string representation
  void output(std::ostream &stream) const;
  
private:
  std::shared_ptr<Vec> values_;    ///< vector of all values, for every components the first value, then the 2nd for all components, etc.
  int nComponents_;    ///< number of components for this field variable, important for interpreting values_
  int componentIndex_; ///< index of the current component for this field variable, starts with 0, important for interpreting values_
 
  std::string name_;    ///< identifier of the component, e.g. 'x'
  std::string exfileBasisFunctionSpecification_;   ///< the basis function specification in the exelem file, e.g. c.Hermite*c.Hermite*c.Hermite
  element_no_t nElements_;      ///< number of elements
  std::shared_ptr<ElementToDofMapping> elementToDofMapping_;   ///< mapping from element-local dof-indices to dof numbers
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping_;   ///< mapping from nodes to dof numbers
  std::shared_ptr<ExfileRepresentation> exfileRepresentation_; ///< indexing for exelem file
};

// output operator
template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const Component<BasisOnMeshType> &rhs);

};  // namespace
#include "field_variable/unstructured/component.tpp"
