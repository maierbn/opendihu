#include "quadrature/newton_cotes.h"

#include <Vc/Vc>
#include <cmath>

namespace Quadrature
{
// 1 NewtonCotes point --> rectangle
template<>
Vc::array<double, 1> NewtonCotes<1>::
samplingPoints()
{
  return Vc::array<double, 1>{
    0.5
  };
}

// 2 NewtonCotes points --> trapezoidal
template<>
Vc::array<double, 2> NewtonCotes<2>::
samplingPoints()
{
  return Vc::array<double, 2>{
    0.,
    1.
  };
}

// 3 NewtonCotes points --> simpson
template<>
Vc::array<double, 3> NewtonCotes<3>::
samplingPoints()
{
  return Vc::array<double, 3>{
    0., 
      0.5, 
      1.
  };
}

// 4 NewtonCotes points --> 3/8
template<>
Vc::array<double, 4> NewtonCotes<4>::
samplingPoints()
{
  return Vc::array<double, 4>{
    0., 
      (1. / 3.), 
      (2. / 3.), 
      1.
  };
}

// 5 NewtonCotes points
template<>
Vc::array<double, 5> NewtonCotes<5>::
samplingPoints()
{
  return Vc::array<double, 5>{
    0.,
      1./4.,
      2./4.,
      3./4.,
      1.
  };
}

// 6 NewtonCotes points
template<>
Vc::array<double, 6> NewtonCotes<6>::
  samplingPoints()
{
  return Vc::array<double, 6>{
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
Vc::array<double, 7> NewtonCotes<7>::
samplingPoints()
{
  return Vc::array<double, 7>{
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
Vc::array<double, 8> NewtonCotes<8>::
samplingPoints()
{
  return Vc::array<double, 8>{
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

// --------------------
// quadrature weights

// 1 NewtonCotes point --> rectangle
template<>
const Vc::array<double, 1> NewtonCotes<1>::
quadratureWeights()
{
  return Vc::array<double, 1>{
    1.0
  };
}

// 2 NewtonCotes points --> trapezoidal
template<>
const Vc::array<double, 2> NewtonCotes<2>::
quadratureWeights()
{
  return Vc::array<double, 2>{
    1./2,
    1./2
  };
}

// 3 NewtonCotes points --> simpson
template<>
const Vc::array<double, 3> NewtonCotes<3>::
quadratureWeights()
{
  return Vc::array<double, 3>{
    1./6.,
    2./3., 
    1./6.
  };
}

// 4 NewtonCotes points --> 3/8
template<>
const Vc::array<double, 4> NewtonCotes<4>::
quadratureWeights()
{
  return Vc::array<double, 4>{
    1./8.,
    3./8.,
    3./8.,
    1./8.
  };
}

// 5 NewtonCotes points
template<>
const Vc::array<double, 5> NewtonCotes<5>::
quadratureWeights()
{
  return Vc::array<double, 5>{
    7./90.,
    16./45.,
    2./15.,
    16./45.,
    7./90.
  };
}

// 6 NewtonCotes points
template<>
const Vc::array<double, 6> NewtonCotes<6>::
  quadratureWeights()
{
  return Vc::array<double, 6>{
    19./288.,
    25./96.,
    25./144.,
    25./144.,
    25./96.,
    19./288.
  };
}

// 7 NewtonCotes points
template<>
const Vc::array<double, 7> NewtonCotes<7>::
quadratureWeights()
{
  return Vc::array<double, 7>{
    41./840.,
    9./35.,
    9./280.,
    34./105.,
    9./280.,
    9./35.,
    41./840.
  };
}

// 8 NewtonCotes points
template<>
const Vc::array<double, 8> NewtonCotes<8>::
quadratureWeights()
{
  return Vc::array<double, 8>{
    751./17280.,
    3577./17280.,
    49./640., 
    2989./17280.,
    2989./17280., 
    49./640., 
    3577./17280., 
    751./17280.
  };
}

}  // namespace
