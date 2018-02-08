#pragma once

namespace Integrator
{
 
/** This class simply holds two integrator types, used for a mixed formulation.
 */
template<typename LowOrderIntegratorType, typename HighOrderIntegratorType>
class Mixed
{
public:
  typedef LowOrderIntegratorType LowOrderIntegrator;
  typedef HighOrderIntegratorType HighOrderIntegrator;
};

};   // namespace