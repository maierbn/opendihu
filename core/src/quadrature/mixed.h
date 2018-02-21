#pragma once

namespace Quadrature
{
 
/** This class simply holds two integrator types, used for a mixed formulation.
 */
template<typename LowOrderQuadratureType, typename HighOrderQuadratureType>
class Mixed
{
public:
  typedef LowOrderQuadratureType LowOrderQuadrature;
  typedef HighOrderQuadratureType HighOrderQuadrature;
};

};   // namespace