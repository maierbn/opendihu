#include "quadrature/gauss.h"

#include <Vc/Vc>

namespace Quadrature
{

template<unsigned int NumberGaussPoints>
constexpr int Gauss<NumberGaussPoints>::numberEvaluations()
{
  return NumberGaussPoints;
}

template<unsigned int NumberGaussPoints>
template<typename ValueType>
ValueType Gauss<NumberGaussPoints>::
computeIntegral(const typename std::array<ValueType,Gauss<NumberGaussPoints>::numberEvaluations()> &evaluations)
{
  return Gauss<NumberGaussPoints>::template computeIntegral<ValueType>(evaluations.begin());
}

// quadrature with 1 Gauss point
template<>
template<typename ValueType>
ValueType Gauss<1>::
computeIntegral(const typename std::array<ValueType,1>::const_iterator evaluationsIter)
{
  ValueType result{};
  result = *evaluationsIter;
  return result;
}

// quadrature with more than 1 Gauss points
template<unsigned int NumberGaussPoints>
template<typename ValueType>
ValueType Gauss<NumberGaussPoints>::
computeIntegral(const typename std::array<ValueType,Gauss<NumberGaussPoints>::numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  const Vc::array<double,NumberGaussPoints> weights = quadratureWeights();
  for (int i = 0; i < NumberGaussPoints; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

}  // namespace
