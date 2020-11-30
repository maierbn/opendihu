#pragma once

#include "mesh/structured_regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"
#include "mesh/composite.h"

namespace Mesh
{

// define the mesh class for the surface of an element, e.g. for 3D elements the 2D surface
template<typename MeshType>
struct SurfaceMesh
{
  typedef MeshType type;
};

template<int D>
struct SurfaceMesh<StructuredRegularFixedOfDimension<D>>
{
  typedef StructuredRegularFixedOfDimension<D-1> type;
};

template<int D>
struct SurfaceMesh<StructuredDeformableOfDimension<D>>
{
  typedef StructuredDeformableOfDimension<D-1> type;
};

template<int D>
struct SurfaceMesh<UnstructuredDeformableOfDimension<D>>
{
  typedef UnstructuredDeformableOfDimension<D-1> type;
};

template<int D>
struct SurfaceMesh<CompositeOfDimension<D>>
{
  typedef StructuredDeformableOfDimension<D-1> type;
};


} // namespace Mesh
