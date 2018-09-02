#include "basis_function/complete_polynomial.h"

namespace BasisFunction
{

// constant functions
// 1D
template<>
double CompletePolynomialOfDimensionAndOrder<1,0>::
phi(int dofIndex, std::array<double,1> xi)
{
  return 1;
}

// 2D
template<>
double CompletePolynomialOfDimensionAndOrder<2,0>::
phi(int dofIndex, std::array<double,2> xi)
{
  return 1;
}

// 3D
template<>
double CompletePolynomialOfDimensionAndOrder<3,0>::
phi(int dofIndex, std::array<double,3> xi)
{
  return 1;
}

// order 0 dphi_dxi
// 1D
template<>
double CompletePolynomialOfDimensionAndOrder<1,0>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,1> xi)
{
  return 0;
}

// 2D
template<>
double CompletePolynomialOfDimensionAndOrder<2,0>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,2> xi)
{
  return 0;
}

// 3D
template<>
double CompletePolynomialOfDimensionAndOrder<3,0>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi)
{
  return 0;
}

// order 1 phi
// 1D
template<>
double CompletePolynomialOfDimensionAndOrder<1,1>::
phi(int dofIndex, std::array<double,1> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return xi[0];
  }
  return 0.0;   // should not be reached
}

// 2D
template<>
double CompletePolynomialOfDimensionAndOrder<2,1>::
phi(int dofIndex, std::array<double,2> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return xi[0];
  case 2:
    return xi[1];
  }
  return 0.0;   // should not be reached
}

// 3D
template<>
double CompletePolynomialOfDimensionAndOrder<3,1>::
phi(int dofIndex, std::array<double,3> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return xi[0];
  case 2:
    return xi[1];
  case 3:
    return xi[2];
  }
  return 0.0;   // should not be reached
}

// order 1 dphi_dxi
// 1D
template<>
double CompletePolynomialOfDimensionAndOrder<1,1>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,1> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 0;
  case 1:
    return 1;
  }
  return 0.0;   // should not be reached
}

// 2D
template<>
double CompletePolynomialOfDimensionAndOrder<2,1>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,2> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 0;
  case 1:
    return (derivativeIdx == 0? 1 : 0);
  case 2:
    return (derivativeIdx == 1? 1 : 0);
  }
  return 0.0;   // should not be reached
}

// 3D
template<>
double CompletePolynomialOfDimensionAndOrder<3,1>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 0;
  case 1:
    return (derivativeIdx == 0? 1 : 0);
  case 2:
    return (derivativeIdx == 1? 1 : 0);
  case 3:
    return (derivativeIdx == 2? 1 : 0);
  }
  return 0.0;   // should not be reached
}

// order 2 phi
// 1D
template<>
double CompletePolynomialOfDimensionAndOrder<1,2>::
phi(int dofIndex, std::array<double,1> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return xi[0];
  case 2:
    return xi[0]*xi[0];
  }
  return 0.0;   // should not be reached
}

// 2D
template<>
double CompletePolynomialOfDimensionAndOrder<2,2>::
phi(int dofIndex, std::array<double,2> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return xi[0];
  case 2:
    return xi[1];
  case 3:
    return xi[0]*xi[0];
  case 4:
    return xi[0]*xi[1];
  case 5:
    return xi[1]*xi[1];
  }
  return 0.0;   // should not be reached
}

// 3D
template<>
double CompletePolynomialOfDimensionAndOrder<3,2>::
phi(int dofIndex, std::array<double,3> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return xi[0];
  case 2:
    return xi[1];
  case 3:
    return xi[2];
  case 4:
    return xi[0]*xi[0];
  case 5:
    return xi[0]*xi[1];
  case 6:
    return xi[0]*xi[2];
  case 7:
    return xi[1]*xi[1];
  case 8:
    return xi[1]*xi[2];
  case 9:
    return xi[2]*xi[2];
  }
  return 0.0;   // should not be reached
}

// order 2 dphi_dxi
// 1D
template<>
double CompletePolynomialOfDimensionAndOrder<1,2>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,1> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 0;
  case 1:
    return 1;
  case 2:
    return 2*xi[0];
  }
  return 0.0;   // should not be reached
}

// 2D
template<>
double CompletePolynomialOfDimensionAndOrder<2,2>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,2> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 0;
  case 1:
    return (derivativeIdx == 0? 1 : 0);
  case 2:
    return (derivativeIdx == 1? 1 : 0);
  case 3:
    return (derivativeIdx == 0? 2*xi[0] : 0);
  case 4:
    return (derivativeIdx == 0? xi[1] : xi[0]);
  case 5:
    return (derivativeIdx == 1? 2*xi[1] : 0);
  }
  return 0.0;   // should not be reached
}

// 3D
template<>
double CompletePolynomialOfDimensionAndOrder<3,2>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi)
{
  switch (dofIndex)
  {
  case 0:
    return 1;
  case 1:
    return (derivativeIdx == 0? 1 : 0);
  case 2:
    return (derivativeIdx == 1? 1 : 0);
  case 3:
    return (derivativeIdx == 2? 1 : 0);
  case 4:
    return (derivativeIdx == 0? 2*xi[0] : 0);
  case 5:
    switch (derivativeIdx)
    {
    case 0:
      return xi[1];
    case 1:
      return xi[0];
    default:
      return 0;
    }
  case 6:
    switch (derivativeIdx)
    {
    case 0:
      return xi[2];
    case 2:
      return xi[0];
    default:
      return 0;
    }
  case 7:
    return (derivativeIdx == 1? 2*xi[1] : 0);
  case 8:
    switch (derivativeIdx)
    {
    case 1:
      return xi[2];
    case 2:
      return xi[1];
    default:
      return 0;
    }
  case 9:
    return (derivativeIdx == 2? 2*xi[2] : 0);
  }
  return 0.0;   // should not be reached
}

};  //namespace
