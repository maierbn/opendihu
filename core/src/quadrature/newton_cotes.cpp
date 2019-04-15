#include "quadrature/newton_cotes.h"

#include <array>
#include <cmath>

namespace Quadrature
{
// 1 NewtonCotes point --> rectangle
template<>
std::array<double, 1> NewtonCotes<1>::
samplingPoints()
{
  return std::array<double, 1>{
    0.5
  };
}

// 2 NewtonCotes points --> trapezoidal
template<>
std::array<double, 2> NewtonCotes<2>::
samplingPoints()
{
  return std::array<double, 2>{
    0.,
    1.
  };
}

// 3 NewtonCotes points --> simpson
template<>
std::array<double, 3> NewtonCotes<3>::
samplingPoints()
{
  return std::array<double, 3>{
    0., 
      0.5, 
      1.
  };
}

// 4 NewtonCotes points --> 3/8
template<>
std::array<double, 4> NewtonCotes<4>::
samplingPoints()
{
  return std::array<double, 4>{
    0., 
      (1. / 3.), 
      (2. / 3.), 
      1.
  };
}

// 5 NewtonCotes points
template<>
std::array<double, 5> NewtonCotes<5>::
samplingPoints()
{
  return std::array<double, 5>{
    0.,
      1./4.,
      2./4.,
      3./4.,
      1.
  };
}

// 6 NewtonCotes points
template<>
std::array<double, 6> NewtonCotes<6>::
  samplingPoints()
{
  return std::array<double, 6>{
    0.,
      1./5.,
      2./5.,
      3./5.,
      4./5.,
      1.,
  };
}

// 7 NewtonCotes points
template<>
std::array<double, 7> NewtonCotes<7>::
samplingPoints()
{
  return std::array<double, 7>{
    0.,
      1./6.,
      2./6.,
      3./6.,
      4./6.,
      5./6.,
      1.,
  };
}

// 8 NewtonCotes points
template<>
std::array<double, 8> NewtonCotes<8>::
samplingPoints()
{
  return std::array<double, 8>{
    0.,
      1. / 7.,
      2. / 7.,
      3. / 7.,
      4. / 7.,
      5. / 7.,
      6. / 7.,
      1.
  };
}

}  // namespace
