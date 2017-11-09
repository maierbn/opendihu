#include "mesh/regular_fixed.h"

#include <array>

namespace Mesh
{

template<>
element_idx_t RegularFixed<1>::nElements()
{
  return nElements_[0];
}

template<>
element_idx_t RegularFixed<2>::nElements()
{
  return nElements_[0]*nElements_[1];
}

template<>
element_idx_t RegularFixed<3>::nElements()
{
  return nElements_[0]*nElements_[1]*nElements_[2];
}


template<>
element_idx_t RegularFixed<1ul>::nDegreesOfFreedom()
{
  return nElements_[0]+1;
}

template<>
element_idx_t RegularFixed<2ul>::nDegreesOfFreedom()
{
  return (nElements_[0]+1) * (nElements_[1]+1);
}

template<>
element_idx_t RegularFixed<3ul>::nDegreesOfFreedom()
{
  return (nElements_[0]+1) * (nElements_[1]+1) * (nElements_[2]+1);
}

};