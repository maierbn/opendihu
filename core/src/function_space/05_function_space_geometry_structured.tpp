#include "function_space/05_function_space_geometry.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType>
std::shared_ptr<typename FunctionSpaceGeometryData<MeshType,BasisFunctionType>::FieldVariableBaseFunctionSpaceType>
FunctionSpaceGeometryData<MeshType,BasisFunctionType>::
fieldVariable(std::string name)
{
  return nullptr;
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceGeometryData<MeshType,BasisFunctionType>::
getDofNoLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, int nodalDofIndex, bool &isOnLocalDomain)
{
  std::array<int,MeshType::dim()> coordinatesLocal = this->meshPartition_->getCoordinatesLocal(coordinatesGlobal, isOnLocalDomain);

  if (isOnLocalDomain)
  {
    node_no_t nodeNoLocal = this->getNodeNo(coordinatesLocal);
    dof_no_t dofNoLocal = nodeNoLocal * FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode() + nodalDofIndex;
    return dofNoLocal;
  }
  return -1;
}

} // namespace
