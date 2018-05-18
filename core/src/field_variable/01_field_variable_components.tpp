#include "field_variable/01_field_variable_components.h"

namespace FieldVariable
{

template<typename BasisOnMeshType,int nComponents>
const std::array<std::string,nComponents> &FieldVariableComponents<BasisOnMeshType,nComponents>::
componentNames() const
{
  return componentNames_;
}

template<typename BasisOnMeshType,int nComponents>
const std::string FieldVariableComponents<BasisOnMeshType,nComponents>::
componentName(int componentNo) const
{
  return componentNames_[componentNo];
}

template<typename BasisOnMeshType, int nComponents>
int FieldVariableComponents<BasisOnMeshType,nComponents>::
getNComponents() const
{
  return nComponents;
}

template<typename BasisOnMeshType, int nComponents>
int FieldVariableComponents<BasisOnMeshType,nComponents>::
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
