
#include "mesh/structured.h"

namespace Mesh
{

template<>
element_no_t Structured<1>::nLocalElements() const
{
  return nElementsPerCoordinateDirectionLocal_[0];
}

template<>
element_no_t Structured<2>::nLocalElements() const
{
  return nElementsPerCoordinateDirectionLocal_[0]*nElementsPerCoordinateDirectionLocal_[1];
}

template<>
element_no_t Structured<3>::nLocalElements() const
{
  return nElementsPerCoordinateDirectionLocal_[0]*nElementsPerCoordinateDirectionLocal_[1]*nElementsPerCoordinateDirectionLocal_[2];
}


}; // namespace