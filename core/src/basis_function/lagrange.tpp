#include "basis_function/lagrange.h"

namespace BasisFunction
{
  
template<int order>
constexpr int LagrangeOfOrder<order>::
nDofsPerNode()
{
  return 1;
}

// linear Lagrange basis
template<>
constexpr int LagrangeOfOrder<1>::
nDofsPerBasis()
{
  return 2;
}

// quadratic Lagrange basis
template<>
constexpr int LagrangeOfOrder<2>::
nDofsPerBasis()
{
  return 3;
}

};  // namespace