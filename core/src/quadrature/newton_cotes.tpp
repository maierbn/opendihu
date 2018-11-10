// all weights are taken from: https://de.wikipedia.org/wiki/Newton-Cotes-Formeln

#include "quadrature/newton_cotes.h"

#include <array>

namespace Quadrature
{

template<unsigned int NumberIntegrationPoints>
constexpr int NewtonCotes<NumberIntegrationPoints>::numberEvaluations()
{
  return NumberIntegrationPoints;
}

template<unsigned int NumberIntegrationPoints>
template<typename ValueType>
ValueType NewtonCotes<NumberIntegrationPoints>::computeIntegral(const typename std::array<ValueType,numberEvaluations()> &evaluations)
{
  return NewtonCotes<NumberIntegrationPoints>::template computeIntegral<ValueType>(evaluations.begin());
}

// quadrature with 1 NewtonCotes point: rectangle rule
template<>
template<typename ValueType>
ValueType NewtonCotes<1>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  result = *evaluationsIter;
  return result;
}
  
// quadrature with 2 NewtonCotes points: trapezoidal rule
template<>
template<typename ValueType>
ValueType NewtonCotes<2>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,2> weights {1./2,1./2};
  for (int i = 0; i < 2; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// quadrature with 3 NewtonCotes points: Simpson's rule
template<>
template<typename ValueType>
ValueType NewtonCotes<3>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,3> weights {1./6., 2./3., 1./6.};
  for (int i = 0; i < 3; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 4 NewtonCotes points: 3/8-rule
template<>
template<typename ValueType>
ValueType NewtonCotes<4>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,4> weights {
    1./8., 3./8., 3./8., 1./8.
  };
  for (int i = 0; i < 4; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 5 NewtonCotes points
template<>
template<typename ValueType>
ValueType NewtonCotes<5>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,5> weights {
    7./90., 16./45., 2./15., 16./45., 7./90.
  };
  for (int i = 0; i < 5; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 6 NewtonCotes points
template<>
template<typename ValueType>
ValueType NewtonCotes<6>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,6> weights {
    19./288., 25./96., 25./144., 25./144., 25./96., 19./288.
  };
  for (int i = 0; i < 6; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 7 NewtonCotes points
template<>
template<typename ValueType>
ValueType NewtonCotes<7>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,7> weights {
    41./840., 9./35., 9./280., 34./105., 9./280., 9./35., 41./840.
  };
  for (int i = 0; i < 7; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}
  
// integrate with 8 NewtonCotes points
template<>
template<typename ValueType>
ValueType NewtonCotes<8>::
computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,8> weights {
    751./17280., 3577./17280., 49./640., 2989./17280., 2989./17280., 49./640., 3577./17280., 751./17280.
  };
  for (int i = 0; i < 8; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}
  
};
