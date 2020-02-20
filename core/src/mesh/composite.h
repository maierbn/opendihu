#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

#include "control/types.h"
#include "mesh/mesh.h"

namespace Mesh
{

/**
 * A composite mesh is a mesh that consists of multiple submeshes of dimension D.
 * Using this it is possible to model a geometry other than a cuboid
 * that behaves the same way as a normal structured mesh.
 * @param D dimension
 */
template<int D>
class CompositeOfDimension : public MeshOfDimension<D>
{
public:
  //! constructor from python settings, one for each submesh
  CompositeOfDimension(PythonConfig specificSettings);

  //! get the total number of elements in the local domain
  element_no_t nElementsLocal() const;

  //! get the total number of elements in the global domain
  global_no_t nElementsGlobal() const;

protected:

  element_no_t nElementsLocal_;     // local number of elements which is the sum of local elements of all submeshes
  global_no_t nElementsGlobal_;     // global number of elements which is the sum of global elements of all submeshes
};

}  // namespace

#include "mesh/composite.tpp"

