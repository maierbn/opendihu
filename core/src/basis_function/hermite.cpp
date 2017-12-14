#include "basis_function/hermite.h"

namespace BasisFunction
{

constexpr int Hermite::nDofsPerNode()
{
  return 2;
}

constexpr int Hermite::nDofsPerBasis()
{
  return 4;
}

constexpr int Hermite::averageNDofsPerElement()
{
  return 2;
}

double Hermite::dphi_dxi(int i, double xi)
{
  if (i==0)
    return -6*xi + 6*xi*xi;
  else if (i==1)
    return 3*xi*xi - 4*xi + 1;
  else if (i==2)
    return -6*xi*xi + 6*xi;
  return 3*xi*xi - 2*xi;
}

double Hermite::phi(int i, double xi)
{
  if (i==0)
    return 1 - 3*xi*xi + 2*xi*xi*xi;
  else if (i==1)
    return xi * (xi-1) * (xi-1);      // xi*xi*xi - 2*xi*xi + xi
  else if (i==2)
    return xi*xi * (3 - 2*xi);        // -2*xi*xi*xi + 3*xi*xi
  return xi*xi * (xi-1);              // xi*xi*xi - xi*xi
}

};  //namespace