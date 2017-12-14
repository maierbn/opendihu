#pragma once

#include <array>

#include <Python.h>

#include "control/types.h"
#include "mesh/mesh.h"

namespace Mesh
{
  
template<unsigned long D>
class Regular : public MeshD<D>
{
public:
  //! constructor
  Regular(std::array<element_idx_t, D> &nElements);
  Regular(PyObject *specificSettings);
  
  //! get number of elements
  element_idx_t nElements(int dimension);
  element_idx_t nElements();
  
  //! get the number of nodes
  element_idx_t nNodes(int dimension);
  element_idx_t nNodes();

protected:
 
  std::array<element_idx_t, D> nElements_;    ///< the number of elements in each coordinate direction
};

};    // namespace

#include "mesh/regular.tpp"