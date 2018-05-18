#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>
#include <petscvec.h>
#include "easylogging++.h"

#include "field_variable/01_field_variable_components.h"

namespace FieldVariable
{

// forward declaration
template<typename BasisOnMeshType, int nComponents> class FieldVariable;

/** Field variable for a structured mesh, i.e. dof and node information are purely implicit.
 *  This is used for RegularFixed and StructuredDeformable meshes.
 */
template<typename BasisOnMeshType, int nComponents_>
class FieldVariableDataStructured :
  public FieldVariableComponents<BasisOnMeshType,nComponents_>
{
public:
  //! inherited constructor
  using FieldVariableComponents<BasisOnMeshType,nComponents_>::FieldVariableComponents;

  //! empty contructor
  FieldVariableDataStructured();

  //! contructor as data copy with a different name (component names are the same)
  FieldVariableDataStructured(FieldVariable<BasisOnMeshType,nComponents_> &rhs, std::string name);

  //! constructor with mesh, name and components
  FieldVariableDataStructured(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames);

  //! destructor
  virtual ~FieldVariableDataStructured();

  //! set all data but the values from a second field variable
  template<typename FieldVariableType>
  void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames);

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

  //! get the number of dofs, i.e. the number of entries per component
  dof_no_t nDofs() const;

  //! if the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate"
  bool isGeometryField() const;



  //! not implemented interface methods

  //! parse current component's exfile representation from file contents
  virtual void parseHeaderFromExelemFile(std::string content){}

  //! parse single element from exelem file
  virtual void parseElementFromExelemFile(std::string content){}

  //! read in values frorm exnode file
  virtual void parseFromExnodeFile(std::string content){}

  //! resize internal representation variable to number of elements
  virtual void setNumberElements(element_no_t nElements){}

  //! reduce memory consumption by removing duplicates in ExfileRepresentations
  virtual void unifyMappings(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode){}

  //! eliminate duplicate elementToDof and exfileRepresentation objects in components of two field variables (this and one other)
  virtual void unifyMappings(std::shared_ptr<FieldVariableBase<BasisOnMeshType>> fieldVariable2){}

  //! initialize PETSc vector with size of total number of dofs for all components of this field variable
  virtual void initializeValuesVector(){}

  //! return the component by index
  virtual std::shared_ptr<Component<BasisOnMeshType>> component(int componentNo) {return nullptr;}   // return empty Component

  //! get the element to dof mapping object
  virtual std::shared_ptr<ElementToDofMapping> elementToDofMapping() const {return nullptr;}

  //! get the node to dof mapping object
  virtual std::shared_ptr<NodeToDofMapping> nodeToDofMapping() const {return nullptr;}

  //! get the number of scale factors
  virtual int getNumberScaleFactors(element_no_t elementGlobalNo) const {return 0;}


protected:

  //! get the number of entries of the internal values_ Vector
  std::size_t nEntries() const;

  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
  std::size_t nEntries_;       ///< number of entries the PETSc vector values_ will have (if it is used). This number of dofs * nComponents

  Vec values_ = PETSC_NULL;          ///< Petsc vector containing the values, the values for the components are stored contiguous, e.g. (comp1val1, comp2val1, comp3val1, comp1val2, comp2val2, ...). Dof ordering proceeds fastest over dofs of a node, then over nodes, node numbering is along whole domain, fastes in x, then in y,z direction.
};

};  // namespace

#include "field_variable/structured/02_field_variable_data_structured.tpp"
