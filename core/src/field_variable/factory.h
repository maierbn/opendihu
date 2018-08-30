#pragma once

#include "field_variable/00_field_variable_base.h"

namespace FieldVariable
{

/** To create a field variable the number of components must be known at compile time.
 *  However, that is not possible when parsing input files. Therefore this class is able to create
 *  field variables with given number of components at runtime.
  */
template<typename FunctionSpaceType>
class Factory
{
public:
  //! create a new empty FieldVariable object with the specified number of components.
  template <typename ...Args>
  static std::shared_ptr<FieldVariableBase<FunctionSpaceType>> makeShared(const int nComponents, Args && ...args);

  //! create a field variable as copy of fieldVariable, with the give components
  template<typename FieldVariableType>
  static std::shared_ptr<FieldVariableBase<FunctionSpaceType>> createFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames);

};

};  // namespace

#include "field_variable/factory.tpp"