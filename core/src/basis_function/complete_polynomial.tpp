#pragma once

#include "basis_function/complete_polynomial.h"

namespace BasisFunction
{
  
template<int D, int order>
constexpr int CompletePolynomialOfDimensionAndOrder<D,order>::
nDofsPerNode()
{
  return 0;
}

// 1D
template<int order>
constexpr int CompletePolynomialNDofs<1,order>::
nDofsPerBasis()
{
  return order+1;
}

// 2D
template<int order>
constexpr int CompletePolynomialNDofs<2,order>::
nDofsPerBasis()
{
  return (order+1)*(order+2)/2;
}

// 3D
template<int order>
constexpr int CompletePolynomialNDofs<3,order>::
nDofsPerBasis()
{
  return (order+1)*(order+2)*(order+3)/6;
}

};  // namespace