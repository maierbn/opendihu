#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "field_variable/00_field_variable_base.h"
#include "field_variable/unstructured/component.h"

namespace FieldVariable
{

/** Base class for a field variable that also stores the component names
 */
template<typename FunctionSpaceType, int nComponentsValue>
class FieldVariableComponents :
  public FieldVariableBaseFunctionSpace<FunctionSpaceType>
{
public:
  //! inherited constructor
  using FieldVariableBaseFunctionSpace<FunctionSpaceType>::FieldVariableBaseFunctionSpace;

  //! get the component names
  const std::array<std::string,nComponentsValue> &componentNames() const;

  //! return the component by index
  virtual std::shared_ptr<Component<FunctionSpaceType,nComponentsValue>> component(int componentNo) = 0;

  //! get the component Name
  const std::string componentName(int componentNo) const;

  //! get the number of components
  static constexpr int nComponents();

  //! get the number of components
  virtual int getNComponents() const;

protected:

  //! get the componentNo that matches the componentName
  int findComponent(std::string componentName);

  std::array<std::string,nComponentsValue> componentNames_;   ///< names of the components, e.g. "x","y","z"
};

}  // namespace
#include "field_variable/01_field_variable_components.tpp"
