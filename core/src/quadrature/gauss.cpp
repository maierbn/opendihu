#include "quadrature/gauss.h"

#include <array>
#include <cmath>

namespace Quadrature
{
  
// 1 Gauss point
template<>
std::array<double, 1> Gauss<1>::
samplingPoints()
{
  return std::array<double, 1>{0.5};
}

// 2 Gauss points
template<>
std::array<double, 2> Gauss<2>::
samplingPoints()
{
  return std::array<double, 2>{
    (-1./sqrt(3.)+1)/2.,
    (+1./sqrt(3.)+1)/2.,
  };
}

// 3 Gauss points
template<>
std::array<double, 3> Gauss<3>::
samplingPoints()
{
  return std::array<double, 3>{
    (-sqrt(3./5)+1)/2.,
    1./2,
    (+sqrt(3./5)+1)/2.,
  };
}

// 4 Gauss points
template<>
std::array<double, 4> Gauss<4>::
samplingPoints()
{
  return std::array<double, 4>{
    (-sqrt(3./7+2./7*sqrt(6./5))+1)/2.,
    (-sqrt(3./7-2./7*sqrt(6./5))+1)/2.,
    (+sqrt(3./7-2./7*sqrt(6./5))+1)/2.,
    (+sqrt(3./7+2./7*sqrt(6./5))+1)/2.
  };
}

// 5 Gauss points
template<>
std::array<double, 5> Gauss<5>::
samplingPoints()
{
  return std::array<double, 5>{
    (-1./3*sqrt(5+2.*sqrt(10./7))+1)/2.,
    (-1./3*sqrt(5-2.*sqrt(10./7))+1)/2.,
    1./2.,
    (+1./3*sqrt(5-2.*sqrt(10./7))+1)/2.,
    (+1./3*sqrt(5+2.*sqrt(10./7))+1)/2.
  };
}

// 7 Gauss points
template<>
std::array<double,6> Gauss<6>::
samplingPoints()
{
  return std::array<double,6>{
    (0.6612093864662645+1)/2.,
    (-0.6612093864662645+1)/2.,
    (-0.2386191860831969+1)/2.,
    (0.2386191860831969+1)/2.,
    (-0.9324695142031521+1)/2.,
    (0.9324695142031521+1)/2.,
  };
}

// 7 Gauss points
template<>
std::array<double,7> Gauss<7>::
samplingPoints()
{
  return std::array<double,7>{
    (0.0000000000000000+1)/2.,
    (0.4058451513773972+1)/2.,
    (-0.4058451513773972+1)/2.,
    (-0.7415311855993945+1)/2.,
    (0.7415311855993945+1)/2.,
    (-0.9491079123427585+1)/2.,
    (0.9491079123427585+1)/2.
  };
}

};