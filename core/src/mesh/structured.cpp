
#include "mesh/structured.h"

namespace Mesh
{
  
template<>
element_no_t Structured<1>::nElements() const
{
  return nElementsPerDimension_[0];
}

template<>
element_no_t Structured<2>::nElements() const
{
  return nElementsPerDimension_[0]*nElementsPerDimension_[1];
}

template<>
element_no_t Structured<3>::nElements() const
{
  return nElementsPerDimension_[0]*nElementsPerDimension_[1]*nElementsPerDimension_[2];
}


}; // namespace