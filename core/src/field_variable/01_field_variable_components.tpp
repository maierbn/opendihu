#include "field_variable/01_field_variable_components.h"

namespace FieldVariable
{

template<typename FunctionSpaceType,int nComponents>
const std::array<std::string,nComponents> &FieldVariableComponents<FunctionSpaceType,nComponents>::
componentNames() const
{
  return componentNames_;
}

template<typename FunctionSpaceType,int nComponents>
const std::string FieldVariableComponents<FunctionSpaceType,nComponents>::
componentName(int componentNo) const
{
  return componentNames_[componentNo];
}

template<typename FunctionSpaceType, int nComponentsValue>
constexpr int FieldVariableComponents<FunctionSpaceType,nComponentsValue>::
nComponents()
{
  return nComponentsValue;
}

template<typename FunctionSpaceType, int nComponentsValue>
int FieldVariableComponents<FunctionSpaceType,nComponentsValue>::
getNComponents() const
{
  return nComponentsValue;
}

template<typename FunctionSpaceType, int nComponents>
int FieldVariableComponents<FunctionSpaceType,nComponents>::
findComponent(std::string componentName)
{
  for (int componentNo = 0; componentNo < componentNames_.size(); componentNo++)
  {
    if (componentNames_[componentNo] == componentName)
      return componentNo;
  }
  return 0;
}

} // namespace
