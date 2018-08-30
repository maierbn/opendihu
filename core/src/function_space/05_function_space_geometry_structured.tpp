#include "function_space/05_function_space_geometry.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType>
std::shared_ptr<typename FunctionSpaceGeometryData<MeshType,BasisFunctionType>::FieldVariableBaseType>
FunctionSpaceGeometryData<MeshType,BasisFunctionType>::
fieldVariable(std::string name)
{
  return nullptr;
}

};  // namespace
