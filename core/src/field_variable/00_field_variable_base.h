#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>
#include "field_variable/interface.h"
#include "basis_on_mesh/basis_on_mesh.h"

namespace FieldVariable
{

// forward declarations
class ElementToNodeMapping;
class NodeToDofMapping;
class ElementToDofMapping;

/** Base class for a field variable that just stores the mesh the field variable is defined on.
 */
template<typename BasisOnMeshType>
class FieldVariableBase :
  public Interface<BasisOnMeshType>
{
public:
  FieldVariableBase();

  //! return the mesh of this field variable
  std::shared_ptr<BasisOnMeshType> mesh();

  //! get the name of the field variable
  std::string name() const;

  //! if the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate"
  bool isGeometryField() const;

  //! set the mesh
  virtual void setMesh(std::shared_ptr<BasisOnMeshType> mesh) {LOG(FATAL) << "this is unused";}
  
  //! get the number of components
  virtual int getNComponents() const = 0;

  //! parse current component's exfile representation from file contents
  virtual void parseHeaderFromExelemFile(std::string content) = 0;

  //! parse single element from exelem file
  virtual void parseElementFromExelemFile(std::string content) = 0;

  //! read in values frorm exnode file
  virtual void parseFromExnodeFile(std::string content) = 0;

  //! resize internal representation variable to number of elements
  virtual void setNumberElements(element_no_t nElements) = 0;

  //! reduce memory consumption by removing duplicates in ExfileRepresentations
  virtual void unifyMappings(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode) = 0;

  //! eliminate duplicate elementToDof and exfileRepresentation objects in components of two field variables (this and one other)
  virtual void unifyMappings(std::shared_ptr<FieldVariableBase<BasisOnMeshType>> fieldVariable2) = 0;

  //! initialize PETSc vector with size of total number of dofs for all components of this field variable
  virtual void initializeValuesVector() = 0;

  //! get the component Name
  virtual const std::string componentName(int componentNo) const = 0;

  //! get the element to dof mapping object
  virtual std::shared_ptr<ElementToDofMapping> elementToDofMapping() const = 0;

  //! get the node to dof mapping object
  virtual std::shared_ptr<NodeToDofMapping> nodeToDofMapping() const = 0;

  //! get the number of scale factors
  virtual int getNumberScaleFactors(element_no_t elementGlobalNo) const = 0;

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  virtual bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2) = 0;

  //! write a exelem file header to a stream, for a particular element
  virtual void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo=-1) = 0;

  //! write a exelem file header to a stream, for a particular node
  virtual void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1) = 0;

  //! for a specific component, get a single value from global dof no.
  virtual double getValue(int componentNo, node_no_t dofGlobalNo) = 0;

  //! set all values to a specific value
  virtual void setValues(double value) = 0;

protected:
 
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information
 
  std::shared_ptr<BasisOnMeshType> mesh_;  ///< the mesh for which the field variable is defined
  std::string name_;     ///< name of the field variable
};

// output operator
template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const FieldVariableBase<BasisOnMeshType> &rhs)
{
  stream << rhs.name();
  return stream;
}

}; // namespace
#include "field_variable/00_field_variable_base.tpp"
