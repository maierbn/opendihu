#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>
#include "field_variable/interface.h"
//#include "function_space/function_space.h"

namespace FieldVariable
{

// forward declarations
class ElementToNodeMapping;
class NodeToDofMapping;
class ElementToDofMapping;

class FieldVariableBase
{};

/** Base class for a field variable that just stores the mesh the field variable is defined on.
 */
template<typename FunctionSpaceType>
class FieldVariableBaseFunctionSpace :
  public FieldVariableBase,
  public Interface<FunctionSpaceType>
{
public:
  FieldVariableBaseFunctionSpace();

  //! return the functionSpace of this field variable
  std::shared_ptr<FunctionSpaceType> functionSpace();

  //! get the name of the field variable
  std::string name() const;

  //! set the name of the field variable
  void setName(std::string name);

  //! if the field has the flag "geometry field", i.e. in the exelem file its type was specified as "coordinate"
  bool isGeometryField() const;

  //! set the field variable to be a "geometry field"
  void setGeometryField(bool isGeometryField);

  //! check if there are NaNs or high values in the current variable, if yes output a warning
  void checkNansInfs(int componentNo = 0) const;

  //! get the number of dofs
  dof_no_t nDofsLocalWithoutGhosts() const;

  //! get the number of global dofs
  dof_no_t nDofsGlobal() const;

  //! get the internal PETSc values vector
  virtual Vec &valuesLocal(int componentNo = 0) = 0;

  //! get the internal PETSc values vector
  virtual Vec &valuesGlobal(int componentNo = 0) = 0;

  //! fill a contiguous vector with all components after each other, "struct of array"-type data layout.
  //! after manipulation of the vector has finished one has to call restoreValuesContiguous
  virtual Vec &getValuesContiguous() = 0;

  //! copy the values back from a contiguous representation where all components are in one vector to the standard internal format of PartitionedPetscVec where there is one local vector with ghosts for each component.
  //! this has to be called
  virtual void restoreValuesContiguous() = 0;

  //! set the functionSpace
  virtual void setFunctionSpace(std::shared_ptr<FunctionSpaceType> functionSpace) {LOG(FATAL) << "this is unused";}
  
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
  virtual void unifyMappings(std::shared_ptr<FieldVariableBaseFunctionSpace<FunctionSpaceType>> fieldVariable2) = 0;

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
  virtual double getValue(int componentNo, node_no_t dofGlobalNo) const = 0;

  //! set all values to a specific value
  virtual void setValues(double value) = 0;

  //! set all entries to 0.0
  virtual void zeroEntries() = 0;

  //! for a specific component, get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  virtual void getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false) const = 0;

protected:
 
  bool isGeometryField_;     ///< if the type of this FieldVariable is a coordinate, i.e. geometric information

  std::shared_ptr<FunctionSpaceType> functionSpace_;  ///< the mesh/function_space for which the field variable is defined
  std::string name_;     ///< name of the field variable
};

// output operator
template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, const FieldVariableBaseFunctionSpace<FunctionSpaceType> &rhs)
{
  stream << rhs.name();
  return stream;
}

}  // namespace
#include "field_variable/00_field_variable_base.tpp"
