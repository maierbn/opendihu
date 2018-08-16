#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/01_field_variable_components.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/unstructured/element_to_node_mapping.h"
#include "field_variable/unstructured/node_to_dof_mapping.h"
#include "field_variable/field_variable_data.h"
#include "field_variable/unstructured/component.h"
#include "partition/partitioned_petsc_vec.h"

namespace FieldVariable
{

/** FieldVariable class for UnstructuredDeformable mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableComponents<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableComponents<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableComponents;

  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;

  //! contructor as data copy with a different name (component names are the same)
  FieldVariableData(FieldVariable<BasisOnMeshType,nComponents> &rhs, std::string name);

  //! constructor with mesh, name and components
  FieldVariableData(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames, dof_no_t nDofsPerComponent);

  //! empty constructor
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

  //! get the number of elements
  element_no_t nLocalElements() const;

  //! get the number of nodes
  node_no_t nNodesLocalWithGhosts() const;

  //! get the number of dofs, i.e. the number of entries per component
  dof_no_t nDofsLocalWithGhosts() const;

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

  //! get the internal values vector
  Vec &values();

  //! get the number of scale factors, TODO: local no
  int getNumberScaleFactors(element_no_t elementGlobalNo) const;

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo=-1);

  //! write a exelem file header to a stream, for a particular node
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2);

  friend class BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>;

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
  void unifyMappings(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents2> &fieldVariable);

  //! eliminate duplicate elementToDof and exfileRepresentation objects in components of two field variables (this and one other)
  void unifyMappings(std::shared_ptr<FieldVariableBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> fieldVariable2);

  //! initialize PETSc vector with size of total number of dofs for all components of this field variable
  void initializeValuesVector();

  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the component by name
  Component<BasisOnMeshType,nComponents> &component(std::string name);

  //! return the component by index
  std::shared_ptr<Component<BasisOnMeshType,nComponents>> component(int componentNo);

  //! return the array of components
  std::array<Component<BasisOnMeshType,nComponents>,nComponents> &component();

  //! multiply dof values with scale factors such that scale factor information is completely contained in dof values
  void eliminateScaleFactors();

  //! if the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate"
  bool isGeometryField() const;

  //! output string representation to stream for debugging
  void output(std::ostream &stream) const;

  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  virtual void finishVectorManipulation() = 0;
  
  //! return the internal partitioned petsc vec
  std::shared_ptr<PartitionedPetscVec<BasisOnMeshType,nComponents>> partitionedPetscVec();
  
protected:

  //! create the element to dof mapping at each component
  void createElementToDofMapping(std::shared_ptr<::FieldVariable::ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode);

  //! initialize components
  void initializeComponents(std::vector<std::string> &componentNames, std::string exfileBasisRepresentation);

  int exfileNo_;    ///< number of the fieldvariable in exelem file (index starts at 1)
  std::size_t nEntries_;       ///< number of entries
  element_no_t nElements_;    ///< number of elements
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::array<Component<BasisOnMeshType,nComponents>,nComponents> component_;   ///< one or multiple components of which this field variable consists of. They correspond to the names in this->componentNames_ (derived from FieldVariableComponents)
  std::shared_ptr<ExfileRepresentation> exfileRepresentation_;       ///< the indexing given in the exelem file, this is the same for all components
  std::shared_ptr<ElementToDofMapping> elementToDofMapping_;       ///< the element to dof mapping of all components, this is the same for all components
  std::shared_ptr<ElementToNodeMapping> elementToNodeMapping_;      ///< mapping from element-local node indices to global node numbers
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping_;       ///< the node to dof mapping of all components, this is the same for all components
  std::shared_ptr<PartitionedPetscVec<BasisOnMeshType,nComponents>> values_;     ///< the vector that contains all values, the entries of all components are interleaved, e.g. (val1comp1, val1comp2, val2comp1, val2comp2, ...)
};

};  // namespace

#include "field_variable/unstructured/02_field_variable_data_unstructured_deformable.tpp"
