
#include "mesh/structured.h"

namespace Mesh
{
  
template<>
element_idx_t Structured<1>::nElements() const
{
  return nElements_[0];
}

template<>
element_idx_t Structured<2>::nElements() const
{
  return nElements_[0]*nElements_[1];
}

template<>
element_idx_t Structured<3>::nElements() const
{
  return nElements_[0]*nElements_[1]*nElements_[2];
}


}; // namespace