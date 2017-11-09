#pragma once

#include <array>

#include "mesh/mesh.h"
#include "control/types.h"

namespace Mesh
{

template<unsigned long D>
class RegularFixed : public MeshD<D>
{
public:
 
  //! construct mesh from python settings
  RegularFixed(PyObject *specificSettings);
  
  //! construct mesh directly from values
  RegularFixed(std::array<element_idx_t, D> nElements, std::array<double, D> physicalExtent);
 
  //! get number of elements
  element_idx_t nElements(int dimension);
  element_idx_t nElements();
  
  //! get mesh width
  double meshWidth(int dimension);
  
  //! get the number of degrees of freedom
  element_idx_t nDegreesOfFreedom();
  
private:
 
  std::array<element_idx_t, D> nElements_;
  std::array<double, D> meshWidth_;
};

}  // namespace

#include "mesh/regular_fixed.tpp"