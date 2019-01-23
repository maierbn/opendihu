#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <fstream>
#include <petscvec.h>
#include "control/types.h"

namespace FieldVariable
{

/** This is the Interface for all field variable types, i.e. it contains the methods that any field variable
 *  independent of the template instanciation provides.
 *  Commented out methods are also part of the interface. They can't be included directly because templates can't be
 *  virtual. Their occurence here is therefore only informative.
 */
template<typename FunctionSpaceType>
struct Interface
{

  //! contructor as data copy with a different name (component names are the same)
  //FieldVariable(const FieldVariable &rhs, std::string name);

  //! constructor with functionSpace, name and components
  //FieldVariable(std::shared_ptr<FunctionSpaceType> functionSpace, std::string name, std::vector<std::string> componentNames);

  //! set all data but the values from a second field variable
  //template<typename FieldVariableType>
  //void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)

  //! copy the values from another field variable of the same type
  //void setValues(FieldVariable &rhs);

  //! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  //template<std::size_t nComponents>
  //void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set values for dofs with one component, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  //template<std::size_t nComponents>
  //void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES);

  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different FunctionSpaceTypes)
  virtual Vec &valuesLocal(int componentNo = 0) = 0;

  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different FunctionSpaceTypes)
  virtual Vec &valuesGlobal(int componentNo = 0) = 0;

  //! get the names of the components
  //virtual std::array<std::string, nComponents> componentNames() const = 0;

  //! get the number of elements
  //virtual element_no_t nElementsLocal() const = 0;

  //! for a specific component, get all values
  //template<int N>
  //void getValues(std::string component, std::array<double,N> &values)

  //! for a specific component, get values from their global dof no.s
  //template<int N>
  //void getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)

  //! get values from their global dof no.s for all components
  //template<int N>
  //void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)

  //! for a specific component, get the values corresponding to all element-local dofs
  //template<int N>
  //void getElementValuesComponent(std::string component, element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values)

  //! get the values corresponding to all element-local dofs for all components
  //void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values)

  //! for a specific component, get a single value from global dof no.
  //virtual double getValue(std::string component, node_no_t dofGlobalNo) = 0;

  //! get a single value from global dof no. for all components
  //template<std::size_t nComponents>
  //std::array<double,nComponents> getValue(node_no_t dofGlobalNo)

  //! write a exelem file header to a stream, for a particular element
  virtual void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo, int fieldVariableNo=-1) = 0;

  //! write a exelem file header to a stream, for a particular node
  virtual void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1) = 0;

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  virtual bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2) = 0;

};

} // namespace FieldVariable
