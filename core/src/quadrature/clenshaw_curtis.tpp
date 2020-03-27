// weights are calculated corresponding to the implementation in ClenshawCurtisweights()
// which is located in examples/quadrature/dominik/src/quadrature.cpp

#include "quadrature/clenshaw_curtis.h"

#include <Vc/Vc>

namespace Quadrature
{

template<unsigned int NumberIntegrationPoints>
constexpr int ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()
{
  return NumberIntegrationPoints;
}

template<unsigned int NumberIntegrationPoints>
template<typename ValueType>
ValueType ClenshawCurtis<NumberIntegrationPoints>::
computeIntegral(const typename std::array<ValueType,ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()> &evaluations)
{
  return ClenshawCurtis<NumberIntegrationPoints>::template computeIntegral<ValueType>(evaluations.begin());
}

// quadrature with 1 ClenshawCurtis point
template<>
template<typename ValueType>
ValueType ClenshawCurtis<1>::
computeIntegral(const typename std::array<ValueType,1>::const_iterator evaluationsIter)
{
  ValueType result{};
  result = *evaluationsIter;
  return result;
}

// quadrature with more than 2 integration points
template<unsigned int NumberIntegrationPoints>
template<typename ValueType>
ValueType ClenshawCurtis<NumberIntegrationPoints>::
computeIntegral(const typename std::array<ValueType,ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  const Vc::array<double,NumberIntegrationPoints> weights = quadratureWeights();
  for (int i = 0; i < NumberIntegrationPoints; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}
  
}  // namespace
