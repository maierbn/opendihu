#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <fstream>

namespace FieldVariable
{
 
/** This is the Interface for all field variable types, i.e. it contains the methods that any field variable
 *  independent of the template instanciation provides.
 *  Commented out methods are also part of the interface. They can't be included directly because templates can't be
 *  virtual. Their occurence here is therefore only informative but can be trusted.
 */
template<typename BasisOnMeshType>
struct Interface 
{
  //! set all data but the values from a second field variable
  //template<typename FieldVariableType>
  //void initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
  
  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
  virtual Vec &values() = 0;
  
  //! get the number of components
  virtual int nComponents() const = 0;
  
  //! get the names of the components
  virtual std::vector<std::string> componentNames() const = 0;
  
  //! get the number of elements
  virtual int nElements() const = 0;
  
  //! for a specific component, get values from their global dof no.s
  //template<int N>
  //void getValuesComponent(std::string component, std::array<int,N> dofGlobalNo, std::array<double,N> &values)
  
  //! get values from their global dof no.s for all components
  //template<int N, int nComponents>
  //void getValues(std::array<int,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
    
  //! for a specific component, get the values corresponding to all element-local dofs
  //template<int N>
  //void getElementValuesComponent(std::string component, element_idx_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
  
  //! get the values corresponding to all element-local dofs for all components
  //template<int nComponents>
  //void getElementValues(element_idx_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
  
  //! for a specific component, get a single value from global dof no.
  virtual double getValue(std::string component, node_idx_t dofGlobalNo) = 0;

  //! get a single value from global dof no. for all components
  //template<int nComponents>
  //std::array<double,nComponents> getValue(node_idx_t dofGlobalNo)
  
  //! write a exelem file header to a stream, for a particular element
  virtual void outputHeaderExelem(std::ostream &file, element_idx_t currentElementGlobalNo, int fieldVariableNo=-1) = 0;
  
  //! write a exelem file header to a stream, for a particular node
  virtual void outputHeaderExnode(std::ostream &file, node_idx_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1) = 0;
  
  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  virtual bool haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2) = 0;
  
};

};  // namespace FieldVariable
