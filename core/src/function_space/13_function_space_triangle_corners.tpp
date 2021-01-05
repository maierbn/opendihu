#include "function_space/13_function_space_triangle_corners.h"

#include "partition/mesh_partition/01_mesh_partition_structured.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType,typename Dummy>
double FunctionSpaceTriangleCorners<MeshType,BasisFunctionType,Dummy>::
phi(int dofIndex, std::array<double,MeshType::dim()> xi, element_no_t elementNoLocal) const
{
  // return normal function value
  return FunctionSpaceFunction<MeshType,BasisFunctionType>::phiHexahedralMesh(dofIndex, xi);
}

template<typename MeshType,typename BasisFunctionType,typename Dummy>
double FunctionSpaceTriangleCorners<MeshType,BasisFunctionType,Dummy>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi, element_no_t elementNoLocal) const
{
  return FunctionSpaceFunction<MeshType,BasisFunctionType>::dphi_dxiHexahedralMesh(dofIndex, derivativeIdx, xi);
}

} // namespace
