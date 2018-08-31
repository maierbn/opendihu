#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/01_field_variable_components.h"
#include "function_space/function_space.h"
#include "function_space/06_function_space_dofs_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/unstructured/element_to_node_mapping.h"
#include "field_variable/unstructured/node_to_dof_mapping.h"
#include "field_variable/field_variable_data.h"
#include "field_variable/unstructured/component.h"
#include "partition/partitioned_petsc_vec.h"

namespace FieldVariable
{

/** FieldVariable class for UnstructuredDeformable mesh
 *  The number of elements, nodes and dofs for unstructured meshes is stored by the geometry field of the mesh and not directly at the mesh like for structured meshes.
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableData<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableComponents<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableComponents<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableComponents;

  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! contructor as data copy with a different name (component names are the same)
  FieldVariableData(FieldVariable<FunctionSpaceType,nComponents> &rhs, std::string name);

  //! contructor as data copy with a different name and different components
  template <int nComponents2>
  FieldVariableData(FieldVariable<FunctionSpaceType,nComponents2> &rhs, std::string name, std::vector<std::string> componentNames);

  //! constructor with functionSpace, name and components
  FieldVariableData(std::shared_ptr<FunctionSpaceType> functionSpace, std::string name, std::vector<std::string> componentNames);

  //! empty constructor, this is needed for parsing exfiles
  FieldVariableData();

  //! destructor
  virtual ~FieldVariableData();

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
  
  //! set the internal mesh
  void setFunctionSpace(std::shared_ptr<FunctionSpaceType> functionSpace);

  //! set the property to be geometry field to this field variable
  void setGeometryField(bool isGeometryField=true);

  //! get the number of entries of the internal values_ Vector
  std::size_t nEntries() const;

  //! get the exfile representation object
  std::shared_ptr<ExfileRepresentation> exfileRepresentation() const;

  //! get the element to dof mapping object
  std::shared_ptr<ElementToDofMapping> elementToDofMapping() const;

  //! get the element to dof mapping object
  std::shared_ptr<ElementToNodeMapping> elementToNodeMapping() const;

  //! get the node to dof mapping object
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping() const;

  //! get the number of elements
  element_no_t nElements() const;
  
  //! get the number of nodes
  node_no_t nNodes() const;
  
  //! get the number of dofs
  dof_no_t nDofs() const;
  
  //! get the internal PETSc values vector
  Vec &valuesLocal(int componentNo = 0);

  //! get the internal PETSc values vector
  Vec &valuesGlobal(int componentNo = 0);

  //! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
  //! after manipulation of the vector has finished one has to call restoreContiguousValuesGlobal
  Vec &getContiguousValuesGlobal();

  //! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
  //! this has to be called
  void restoreContiguousValuesGlobal();

  //! get the number of scale factors
  int getNumberScaleFactors(element_no_t elementGlobalNo) const;

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo=-1);

  //! write a exelem file header to a stream, for a particular node
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2);

  //! resize internal representation variable to number of elements
  void setNumberElements(element_no_t nElements);

  //! parse current component's exfile representation from file contents
  void parseHeaderFromExelemFile(std::string content);

  //! parse single element from exelem file
  void parseElementFromExelemFile(std::string content);

  //! read in values frorm exnode file
  void parseFromExnodeFile(std::string content);

  //! reduce memory consumption by removing duplicates in ExfileRepresentations
  void unifyMappings(std::shared_ptr<::FieldVariable::ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode);

  //! eliminate duplicate elementToDof and exfileRepresentation objects in components of two field variables (this and one other)
  template<int nComponents2>
  void unifyMappings(FieldVariable<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents2> &fieldVariable);

  //! eliminate duplicate elementToDof and exfileRepresentation objects in components of two field variables (this and one other)
  void unifyMappings(std::shared_ptr<FieldVariableBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> fieldVariable2);

  //! initialize PETSc vector with size of total number of dofs for all components of this field variable
  void initializeValuesVector();

  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the component by name
  Component<FunctionSpaceType,nComponents> &component(std::string name);

  //! return the component by index
  std::shared_ptr<Component<FunctionSpaceType,nComponents>> component(int componentNo);

  //! return the array of components
  std::array<Component<FunctionSpaceType,nComponents>,nComponents> &component();

  //! multiply dof values with scale factors such that scale factor information is completely contained in dof values
  void eliminateScaleFactors();

  //! output string representation to stream for debugging
  void output(std::ostream &stream) const;

  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  virtual void finishVectorManipulation() = 0;
  
  //! return the internal partitioned petsc vec
  std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents>> partitionedPetscVec();
  
protected:

  //! create the element to dof mapping at each component
  void createElementToDofMapping(std::shared_ptr<::FieldVariable::ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode);

  //! initialize components
  void initializeComponents(std::vector<std::string> &componentNames, std::string exfileBasisRepresentation);

  int exfileNo_;    ///< number of the fieldvariable in exelem file (index starts at 1)
  element_no_t nElements_;    ///< number of elements
  std::array<Component<FunctionSpaceType,nComponents>,nComponents> component_;   ///< one or multiple components of which this field variable consists of. They correspond to the names in this->componentNames_ (derived from FieldVariableComponents)
  std::shared_ptr<ExfileRepresentation> exfileRepresentation_;       ///< the indexing given in the exelem file, this is the same for all components
  std::shared_ptr<ElementToDofMapping> elementToDofMapping_;       ///< the element to dof mapping of all components, this is the same for all components
  std::shared_ptr<ElementToNodeMapping> elementToNodeMapping_;      ///< mapping from element-local node indices to global node numbers
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping_;       ///< the node to dof mapping of all components, this is the same for all components
  std::shared_ptr<PartitionedPetscVec<FunctionSpaceType,nComponents>> values_;     ///< the vector that contains all values, the entries of all components are interleaved, e.g. (val1comp1, val1comp2, val2comp1, val2comp2, ...)
};

};  // namespace

#include "field_variable/unstructured/02_field_variable_data_unstructured_deformable.tpp"
