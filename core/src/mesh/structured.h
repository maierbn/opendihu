#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

#include "control/types.h"
#include "mesh/mesh.h"

namespace Mesh
{

/**
 * A structured mesh, i.e. it has a fixed number of elements in x,y and z direction.
 * This mesh type knows its number of elements.
 */
template<int D>
class Structured : public MeshOfDimension<D>
{
public:
 // TODO: constructor not needed? remove
  //! constructor from number of elements in coordinate directions
  //Structured(std::array<element_no_t, D> &nElements);

  //! constructor from python settings
  Structured(PyObject *specificSettings);

  //! get number of elements in a given coordinate direction for the local domain
  element_no_t nElementsPerCoordinateDirectionLocal(int dimension) const;

  //! get the array with all numbers of elements per coordinate direction for the local domain
  std::array<element_no_t, D> nElementsPerCoordinateDirectionLocal() const;

  //! get number of elements in a given coordinate direction for the global domain
  global_no_t nElementsPerCoordinateDirectionGlobal(int dimension) const;

  //! get the array with all numbers of elements per coordinate direction of the global domain
  std::array<global_no_t, D> nElementsPerCoordinateDirectionGlobal() const;

  //! get the total number of elements in the local domain
  element_no_t nElementsLocal() const;
  
  //! get the total number of elements in the global domain
  global_no_t nElementsGlobal() const;

protected:
  std::array<element_no_t, D> nElementsPerCoordinateDirectionLocal_;    ///< the number of stored elements in each coordinate direction (the locally computed portion)
  std::array<global_no_t, D> nElementsPerCoordinateDirectionGlobal_;    ///< the global number of stored elements in each coordinate direction. This variable is set when reading from config the first time.
  std::array<int, D> nRanks_;    ///< when the number of elements is specified locally the grid of ranks of the partitioning
};

};    // namespace

#include "mesh/structured.tpp"
