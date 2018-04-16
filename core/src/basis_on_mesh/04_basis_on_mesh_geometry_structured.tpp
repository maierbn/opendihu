#include "basis_on_mesh/04_basis_on_mesh_geometry.h"

namespace BasisOnMesh
{
  
template<typename MeshType,typename BasisFunctionType>
std::shared_ptr<typename BasisOnMeshGeometryData<MeshType,BasisFunctionType>::FieldVariableBaseType> 
BasisOnMeshGeometryData<MeshType,BasisFunctionType>::
fieldVariable(std::string name)
{
  return nullptr;
}

};  // namespace
