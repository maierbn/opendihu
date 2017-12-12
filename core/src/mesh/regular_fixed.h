#pragma once

#include <array>

#include "mesh/mesh.h"
#include "mesh/regular.h"
#include "control/types.h"

namespace Mesh
{

template<unsigned long D>
class RegularFixed : public Regular<D>
{
public:
 
  //! construct mesh from python settings
  RegularFixed(PyObject *specificSettings);
  
  //! construct mesh directly from values
  RegularFixed(std::array<element_idx_t, D> nElements, std::array<double, D> physicalExtent);
 
  //! get mesh width
  double meshWidth(int dimension);
  
private:
 
  std::array<double, D> meshWidth_;
};

}  // namespace

#include "mesh/regular_fixed.tpp"