#pragma once

#include "mesh/structured_regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"
#include "mesh/composite.h"

namespace Mesh
{

// structured meshes with given dimension
template<int D,typename MeshType>
using isStructuredWithDim = std::enable_if_t<
  std::is_same<MeshType, StructuredRegularFixedOfDimension<D>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<D>>::value
  ,
  MeshType
>;

// structured meshes
template<typename MeshType>
using isStructured = std::enable_if_t<
  std::is_same<MeshType, StructuredRegularFixedOfDimension<1>>::value
  || std::is_same<MeshType, StructuredRegularFixedOfDimension<2>>::value
  || std::is_same<MeshType, StructuredRegularFixedOfDimension<3>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<1>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<2>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<3>>::value
  ,
  MeshType
>;

// structured and composite meshes
template<typename MeshType>
using isStructuredOrComposite = std::enable_if_t<
  std::is_same<MeshType, StructuredRegularFixedOfDimension<1>>::value
  || std::is_same<MeshType, StructuredRegularFixedOfDimension<2>>::value
  || std::is_same<MeshType, StructuredRegularFixedOfDimension<3>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<1>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<2>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<3>>::value
  || std::is_same<MeshType, CompositeOfDimension<1>>::value
  || std::is_same<MeshType, CompositeOfDimension<2>>::value
  || std::is_same<MeshType, CompositeOfDimension<3>>::value
  ,
  MeshType
>;

// deformable meshes
template<typename MeshType>
using isDeformable = std::enable_if_t<
  std::is_same<MeshType, StructuredDeformableOfDimension<1>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<2>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<3>>::value
  || std::is_same<MeshType, UnstructuredDeformableOfDimension<1>>::value
  || std::is_same<MeshType, UnstructuredDeformableOfDimension<2>>::value
  || std::is_same<MeshType, UnstructuredDeformableOfDimension<3>>::value
  ,
  MeshType
>;

// deformable meshes with given dimension
template<int D,typename MeshType>
using isDeformableWithDim = std::enable_if_t<
  std::is_same<MeshType, StructuredDeformableOfDimension<D>>::value
  || std::is_same<MeshType, UnstructuredDeformableOfDimension<D>>::value
  ,
  MeshType
>;

template<int D,typename MeshType>
using isDim = std::enable_if_t<
  std::is_same<MeshType, StructuredRegularFixedOfDimension<D>>::value
  || std::is_same<MeshType, StructuredDeformableOfDimension<D>>::value
  || std::is_same<MeshType, UnstructuredDeformableOfDimension<D>>::value
  || std::is_same<MeshType, CompositeOfDimension<D>>::value
  ,
  MeshType
>;

template<int D,typename MeshType>
using isNotDim = std::enable_if_t<
  !std::is_same<MeshType, StructuredRegularFixedOfDimension<D>>::value
  && !std::is_same<MeshType, StructuredDeformableOfDimension<D>>::value
  && !std::is_same<MeshType, UnstructuredDeformableOfDimension<D>>::value
  ,
  MeshType
>;

// if field variable has a mesh of type composite
template<typename FieldVariableTypeSharedPtr>
using isComposite = std::is_same<
  typename FieldVariableTypeSharedPtr::element_type::FunctionSpace::Mesh,
  CompositeOfDimension<FieldVariableTypeSharedPtr::element_type::FunctionSpace::Mesh::dim()>
>;

} // namespace Mesh
