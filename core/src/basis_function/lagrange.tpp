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

//! return a basis order string as used in python files and callbacks, e.g. "2"
template<int order>
constexpr int LagrangeOfOrder<order>::
getBasisOrder()
{
  return order;
}

//! return a basis function type string as used in python files and callbacks, e.g. "Lagrange"
template<int order>
std::string LagrangeOfOrder<order>::
getBasisFunctionString()
{
  return std::string("Lagrange");
}

} // namespace
