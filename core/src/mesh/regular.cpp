
#include "mesh/regular.h"

namespace Mesh
{
  
template<>
element_idx_t Regular<1>::nElements()
{
  return nElements_[0];
}

template<>
element_idx_t Regular<2>::nElements()
{
  return nElements_[0]*nElements_[1];
}

template<>
element_idx_t Regular<3>::nElements()
{
  return nElements_[0]*nElements_[1]*nElements_[2];
}


template<>
element_idx_t Regular<1ul>::nNodes()
{
  return nElements_[0]+1;
}

template<>
element_idx_t Regular<2ul>::nNodes()
{
  return (nElements_[0]+1) * (nElements_[1]+1);
}

template<>
element_idx_t Regular<3ul>::nNodes()
{
  return (nElements_[0]+1) * (nElements_[1]+1) * (nElements_[2]+1);
}

}; // namespace