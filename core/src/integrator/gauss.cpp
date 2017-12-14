#include "integrator/gauss.h"

#include <array>
#include <cmath>

namespace Integrator
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

// integrate with 1 Gauss point
template<>
double Gauss<1>::integrate(double *evaluations)
{
  return evaluations[0];
}

// integrate with 2 Gauss points
template<>
double Gauss<2>::integrate(double *evaluations)
{
  std::array<double,2> weights {1./2,1./2};
  double result = 0;
  for (int i=0; i<2; i++)
  {
    result += weights[i]*evaluations[i];
  }
  return result;
}

// integrate with 3 Gauss points
template<>
double Gauss<3>::integrate(double *evaluations)
{
  std::array<double,3> weights {5./18, 4./9., 5./18};
  double result = 0;
  for (int i=0; i<3; i++)
  {
    result += weights[i]*evaluations[i];
  }
  return result;
}

// integrate with 4 Gauss points
template<>
double Gauss<4>::integrate(double *evaluations)
{
  std::array<double,4> weights {
    (18.-sqrt(30.))/36./2.,
    (18.+sqrt(30.))/36./2.,
    (18.+sqrt(30.))/36./2.,
    (18.-sqrt(30.))/36./2.
  };
  double result = 0;
  for (int i=0; i<4; i++)
  {
    result += weights[i]*evaluations[i];
  }
  return result;
}

// integrate with 5 Gauss points
template<>
double Gauss<5>::integrate(double *evaluations)
{
  std::array<double,5> weights {
    (322.-13.*sqrt(70.))/900./2.,
    (322.+13.*sqrt(70.))/900./2.,
    128./225./2.,
    (322.+13.*sqrt(70.))/900./2.,
    (322.-13.*sqrt(70.))/900./2.
  };
  double result = 0;
  for (int i=0; i<5; i++)
  {
    result += weights[i]*evaluations[i];
  }
  return result;
}

// integrate with 6 Gauss points
template<>
double Gauss<6>::integrate(double *evaluations)
{
  std::array<double,6> weights {
    0.3607615730481386/2.,
    0.3607615730481386/2.,
    0.4679139345726910/2.,
    0.4679139345726910/2.,
    0.1713244923791704/2.,
    0.1713244923791704/2.,
  };
  double result = 0;
  for (int i=0; i<6; i++)
  {
    result += weights[i]*evaluations[i];
  }
  return result;
}

// integrate with 7 Gauss points
template<>
double Gauss<7>::integrate(double *evaluations)
{
  std::array<double,7> weights {
    0.4179591836734694/2.,
    0.3818300505051189/2.,
    0.3818300505051189/2.,
    0.2797053914892766/2.,
    0.2797053914892766/2.,
    0.1294849661688697/2.,
    0.1294849661688697/2.,
  };
  double result = 0;
  for (int i=0; i<7; i++)
  {
    result += weights[i]*evaluations[i];
  }
  return result;
}

};