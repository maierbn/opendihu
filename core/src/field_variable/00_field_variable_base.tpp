#include "field_variable/00_field_variable_base.h"

namespace FieldVariable
{

template<typename FunctionSpaceType>
FieldVariableBase<FunctionSpaceType>::
FieldVariableBase() : functionSpace_(nullptr)
{
}

template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> FieldVariableBase<FunctionSpaceType>::
functionSpace()
{
  return functionSpace_;
}

template<typename FunctionSpaceType>
std::string FieldVariableBase<FunctionSpaceType>::
name() const
{
  return this->name_;
}

template<typename FunctionSpaceType>
bool FieldVariableBase<FunctionSpaceType>::
isGeometryField() const
{
  return this->isGeometryField_;
}

} // namespace
