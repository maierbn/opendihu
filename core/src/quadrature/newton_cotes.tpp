// all weights are taken from: https://de.wikipedia.org/wiki/Newton-Cotes-Formeln

#include "quadrature/newton_cotes.h"

#include <Vc/Vc>

namespace Quadrature
{

template<unsigned int NumberIntegrationPoints>
constexpr int NewtonCotes<NumberIntegrationPoints>::numberEvaluations()
{
  return NumberIntegrationPoints;
}

template<unsigned int NumberIntegrationPoints>
template<typename ValueType>
ValueType NewtonCotes<NumberIntegrationPoints>::
computeIntegral(const typename Vc::array<ValueType,NewtonCotes<NumberIntegrationPoints>::numberEvaluations()> &evaluations)
{
  return NewtonCotes<NumberIntegrationPoints>::template computeIntegral<ValueType>(evaluations.begin());
}

// quadrature with 1 NewtonCotes point: rectangle rule
template<>
template<typename ValueType>
ValueType NewtonCotes<1>::
computeIntegral(const typename Vc::array<ValueType,1>::const_iterator evaluationsIter)
{
  ValueType result{};
  result = *evaluationsIter;
  return result;
}

// quadrature with more than 2 integration points
template<unsigned int NumberIntegrationPoints>
template<typename ValueType>
ValueType NewtonCotes<NumberIntegrationPoints>::
computeIntegral(const typename Vc::array<ValueType,NewtonCotes<NumberIntegrationPoints>::numberEvaluations()>::const_iterator evaluationsIter)
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
