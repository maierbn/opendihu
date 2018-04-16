#include "field_variable/factory.h"

#include "easylogging++.h"

namespace FieldVariable
{

template<typename BasisOnMeshType>
template <typename ...Args>
std::shared_ptr<FieldVariableBase<BasisOnMeshType>> Factory<BasisOnMeshType>::
makeShared(const int nComponents, Args && ...args)
{
  switch (nComponents)
  {
    case 1:
      return std::make_shared<FieldVariable<BasisOnMeshType,1>>(std::forward<Args>(args)...);
    case 2:
      return std::make_shared<FieldVariable<BasisOnMeshType,2>>(std::forward<Args>(args)...);
    case 3:
      return std::make_shared<FieldVariable<BasisOnMeshType,3>>(std::forward<Args>(args)...);
    case 4:
      return std::make_shared<FieldVariable<BasisOnMeshType,4>>(std::forward<Args>(args)...);
    default:
      LOG(ERROR) << "Could not create field variable with nComponents=" << nComponents;
      break;
  }
  return nullptr;
};

template<typename BasisOnMeshType>
template<typename FieldVariableType>
std::shared_ptr<FieldVariableBase<BasisOnMeshType>> Factory<BasisOnMeshType>::
createFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
{
  const int nComponents = componentNames.size();
  if (nComponents == 1)
  {
    std::shared_ptr<FieldVariable<BasisOnMeshType,1>> result = std::make_shared<FieldVariable<BasisOnMeshType,1>>();
    result->initializeFromFieldVariable(fieldVariable, name, componentNames);
    return result;
  }
  else if (nComponents == 2)
  {
    std::shared_ptr<FieldVariable<BasisOnMeshType,2>> result = std::make_shared<FieldVariable<BasisOnMeshType,2>>();
    result->initializeFromFieldVariable(fieldVariable, name, componentNames);
    return result;
  }
  else if (nComponents == 3)
  {
    std::shared_ptr<FieldVariable<BasisOnMeshType,3>> result = std::make_shared<FieldVariable<BasisOnMeshType,3>>();
    result->initializeFromFieldVariable(fieldVariable, name, componentNames);
    return result;
  }
  else if (nComponents == 4)
  {
    std::shared_ptr<FieldVariable<BasisOnMeshType,4>> result = std::make_shared<FieldVariable<BasisOnMeshType,4>>();
    result->initializeFromFieldVariable(fieldVariable, name, componentNames);
    return result;
  }
  else
  {
    LOG(ERROR) << "Could not create field variable with nComponents=" << nComponents;
  }
  
  return nullptr;
}
  
};  // namespace
