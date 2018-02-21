#include "integrator/gauss.h"

#include <array>

namespace Integrator
{

template<unsigned int NumberGaussPoints>
constexpr int Gauss<NumberGaussPoints>::numberEvaluations()
{
  return NumberGaussPoints;
}

template<unsigned int NumberGaussPoints>
template<typename ValueType>
ValueType Gauss<NumberGaussPoints>::integrate(const typename std::array<ValueType,numberEvaluations()> &evaluations)
{
  return Gauss<NumberGaussPoints>::template integrate<ValueType>(evaluations.begin());
}

// quadrature with 1 Gauss point
template<>
template<typename ValueType>
ValueType Gauss<1>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  result = *evaluationsIter;
  return result;
}
  
// quadrature with 2 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<2>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  std::array<double,2> weights {1./2,1./2};
  for (int i = 0; i < 2; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// quadrature with 3 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<3>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  std::array<double,3> weights {5./18, 4./9., 5./18};
  for (int i = 0; i < 3; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 4 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<4>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  std::array<double,4> weights {
    (18.-sqrt(30.))/36./2.,
    (18.+sqrt(30.))/36./2.,
    (18.+sqrt(30.))/36./2.,
    (18.-sqrt(30.))/36./2.
  };
  for (int i = 0; i < 4; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 5 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<5>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  std::array<double,5> weights {
    (322.-13.*sqrt(70.))/900./2.,
    (322.+13.*sqrt(70.))/900./2.,
    128./225./2.,
    (322.+13.*sqrt(70.))/900./2.,
    (322.-13.*sqrt(70.))/900./2.
  };
  for (int i = 0; i < 5; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 6 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<6>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  std::array<double,6> weights {
    0.3607615730481386/2.,
    0.3607615730481386/2.,
    0.4679139345726910/2.,
    0.4679139345726910/2.,
    0.1713244923791704/2.,
    0.1713244923791704/2.,
  };
  for (int i = 0; i < 6; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 7 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<7>::
integrate(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result;
  std::array<double,7> weights {
    0.4179591836734694/2.,
    0.3818300505051189/2.,
    0.3818300505051189/2.,
    0.2797053914892766/2.,
    0.2797053914892766/2.,
    0.1294849661688697/2.,
    0.1294849661688697/2.,
  };
  for (int i = 0; i < 7; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}
  
};