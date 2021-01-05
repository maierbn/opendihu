#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <petscmat.h>

#include "control/types.h"
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "mesh/deformable.h"

namespace Mesh
{

/**
 * A structured mesh, i.e. fixed number of elements in each coordinate direction.
 * The shape of the elements/the node positions is arbitrary.
 * The node positions can be changed by computation, e.g. for computational mechanics.
 */
template<int D>
class StructuredDeformableOfDimension :
  public Structured<D>, public Deformable
{
public:
  //! constructor from python settings
  StructuredDeformableOfDimension(PythonConfig specificSettings);

  //! if the mesh has triangles at its corners in each x-y plane, only possible if 3D
  bool hasTriangleCorners();

protected:

  bool hasTriangleCorners_;    //< if the corners in each x-y plane are triangles instead of quadrilaterals, only in 3D and used for the muscle geometry
};

}  // namespace

#include "mesh/structured_deformable.tpp"
