#pragma once

#include "mesh/regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"

namespace Mesh
{

// structured meshes
template<int D,typename MeshType>
using isStructuredWithDim = std::enable_if_t<
  std::is_same<MeshType, RegularFixed<D>>::value
  || std::is_same<MeshType, StructuredDeformable<D>>::value
  ,
  MeshType
>;

// deformable meshes
template<typename MeshType>
using isDeformable = std::enable_if_t<
  std::is_same<MeshType, StructuredDeformable<1>>::value
  || std::is_same<MeshType, StructuredDeformable<2>>::value
  || std::is_same<MeshType, StructuredDeformable<3>>::value
  || std::is_same<MeshType, UnstructuredDeformable<1>>::value
  || std::is_same<MeshType, UnstructuredDeformable<2>>::value
  || std::is_same<MeshType, UnstructuredDeformable<3>>::value
  ,
  MeshType
>;


template<int D,typename MeshType>
using isDim = std::enable_if_t<
  std::is_same<MeshType, RegularFixed<D>>::value
  || std::is_same<MeshType, StructuredDeformable<D>>::value
  || std::is_same<MeshType, UnstructuredDeformable<D>>::value
  ,
  MeshType
>;

} // namespace Mesh
