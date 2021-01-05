#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "mesh/type_traits.h"
#include "function_space/13_function_space_triangle_corners.h"
#include "mesh/mesh.h"
#include "mesh/surface_mesh.h"
#include "basis_function/lagrange.h"
#include "function_space/function_space_generic.h"
#include <Vc/Vc>

namespace FunctionSpace
{

/** FunctionSpace derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace : public FunctionSpaceTriangleCorners<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceTriangleCorners<MeshType,BasisFunctionType>::FunctionSpaceTriangleCorners;

  typedef MeshType Mesh;
  typedef BasisFunctionType BasisFunction;
  typedef FunctionSpace<MeshType,BasisFunctionType> HighOrderFunctionSpace;
  typedef typename ::Mesh::SurfaceMesh<MeshType>::type SurfaceMesh;

  //! get a description of the function space, with mesh name and type
  std::string getDescription() const;
};

}  // namespace

#include "function_space/function_space.tpp"
