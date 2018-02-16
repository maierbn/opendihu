#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "control/types.h"

namespace Mesh
{

/** 
 * A regular mesh with give number of elements in each coordinate direction. The elements
 * have a all the same length/quadratic/cubic size (=mesh width) per coordinate direction.
 * This mesh cannot deform, i.e. it cannot be used for structural mechanics with deformations.
 */
template<int D>
class RegularFixed : public Structured<D>  // StructuredRegularFixed structured_regular_fixed.h
{
public:
  //! constructor from base class
  using Structured<D>::Structured;
 
  //! construct mesh from python settings
  //RegularFixed(PyObject *specificSettings);
  
  //! construct mesh directly from values
  //RegularFixed(std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent);
 
  //! get mesh width of the given coordinate direction
  virtual double meshWidth(int dimension) const = 0;  // defined in field_variable_regular_fixed.tpp
  
private:
 
};

}  // namespace

#include "mesh/regular_fixed.tpp"
