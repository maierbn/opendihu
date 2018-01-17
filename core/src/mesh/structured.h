#pragma once

#include <array>

#include <Python.h>

#include "control/types.h"
#include "mesh/mesh.h"

namespace Mesh
{
  
/**
 * A structured mesh, i.e. it has a fixed number of elements in x,y and z direction.
 * This mesh type knows its number of elements.
 */
template<int D>
class Structured : public MeshD<D>
{
public:
  //! constructor
  Structured(std::array<element_idx_t, D> &nElements);
  Structured(PyObject *specificSettings);
  
  //! get number of elements
  element_idx_t nElements(int dimension) const;
  element_idx_t nElements() const;

protected:
 
  std::array<element_idx_t, D> nElements_;    ///< the number of elements in each coordinate direction
};

};    // namespace

#include "mesh/structured.tpp"