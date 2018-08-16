
#include "mesh/structured.h"

namespace Mesh
{

template<>
element_no_t Structured<1>::nElementsLocal() const
{
  return nElementsPerCoordinateDirectionLocal_[0];
}

template<>
element_no_t Structured<2>::nElementsLocal() const
{
  return nElementsPerCoordinateDirectionLocal_[0]*nElementsPerCoordinateDirectionLocal_[1];
}

template<>
element_no_t Structured<3>::nElementsLocal() const
{
  return nElementsPerCoordinateDirectionLocal_[0]*nElementsPerCoordinateDirectionLocal_[1]*nElementsPerCoordinateDirectionLocal_[2];
}

template<>
global_no_t Structured<1>::nElementsGlobal() const
{
  return nElementsPerCoordinateDirectionGlobal_[0];
}

template<>
global_no_t Structured<2>::nElementsGlobal() const
{
  return nElementsPerCoordinateDirectionGlobal_[0]*nElementsPerCoordinateDirectionGlobal_[1];
}

template<>
global_no_t Structured<3>::nElementsGlobal() const
{
  return nElementsPerCoordinateDirectionGlobal_[0]*nElementsPerCoordinateDirectionGlobal_[1]*nElementsPerCoordinateDirectionGlobal_[2];
}


}; // namespace