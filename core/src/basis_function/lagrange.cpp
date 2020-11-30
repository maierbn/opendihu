#include "basis_function/lagrange.h"

namespace BasisFunction
{

// linear
template<>
double LagrangeOfOrder<1>::
dphi_dxi(int i, double xi)
{
  if (i==0)
    return -1;
  return 1;
}

template<>
double LagrangeOfOrder<1>::
phi(int i, double xi)
{
  if (i==0)
    return 1-xi;
  return xi;
}

// quadratic
template<>
double LagrangeOfOrder<2>::
dphi_dxi(int i, double xi)
{
  if (i==0)
    return 4*xi - 3;
  else if (i==1)
    return -8*xi + 4;
  return 4*xi - 1;
}

template<>
double LagrangeOfOrder<2>::
phi(int i, double xi)
{
  // phi(xi) = a*xi^2 + b*xi + c
  // phi0(xi) = 2*xi^2 -3*xi + 1
  // phi1(xi) = -4*xi^2 + 4*xi
  // phi2(xi) = 2*xi^2 - xi

  if (i==0)
    return (2*xi - 1) * (xi-1);    // 2*xi*xi - 3*xi + 1
  else if (i==1)
    return 4*(xi - xi*xi);
  return 2*xi*xi - xi;
}

}  // namespace
