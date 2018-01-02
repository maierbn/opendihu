#pragma once

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
template<unsigned long D>
class RegularFixed : public Structured<D>
{
public:
 
  //! construct mesh from python settings
  RegularFixed(PyObject *specificSettings);
  
  //! construct mesh directly from values
  RegularFixed(std::array<element_idx_t, D> nElements, std::array<double, D> physicalExtent);
 
  //! get mesh width
  double meshWidth(int dimension) const;
  
private:
 
  std::array<double, D> meshWidth_;
};

}  // namespace

#include "mesh/regular_fixed.tpp"