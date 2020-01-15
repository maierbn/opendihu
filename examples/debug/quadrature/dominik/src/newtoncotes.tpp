// all weights are taken from: https://de.wikipedia.org/wiki/Newton-Cotes-Formeln

#include "newtoncotes.h"

#include <array>

namespace Quadrature
{

template<unsigned int NumberNCPoints>
constexpr int NC<NumberNCPoints>::numberEvaluationsD()
{
  return NumberNCPoints;
}

template<unsigned int NumberNCPoints>
template<typename ValueType>
//ValueType NC<NumberNCPoints>::computeIntegral(const typename std::array<ValueType,numberEvaluations()> &evaluations)
ValueType NC<NumberNCPoints>::computeIntegral(const typename std::array<ValueType,numberEvaluations> &evaluations)
{
  return NC<NumberNCPoints>::template computeIntegral<ValueType>(evaluations.begin());
}

// quadrature with 1 NC point
template<>
template<typename ValueType>
ValueType NC<1>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
{
  ValueType result{};
  result = *evaluationsIter;
  return result;
}
  
// quadrature with 2 NC points
template<>
template<typename ValueType>
ValueType NC<2>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,2> weights {1./2,1./2};
  for (int i = 0; i < 2; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// quadrature with 3 NC points
template<>
template<typename ValueType>
ValueType NC<3>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,3> weights {1./6., 2./3., 1./6.};
  for (int i = 0; i < 3; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 4 NC points
template<>
template<typename ValueType>
ValueType NC<4>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
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

// integrate with 5 NC points
template<>
template<typename ValueType>
ValueType NC<5>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
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

// integrate with 6 NC points
template<>
template<typename ValueType>
ValueType NC<6>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
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

// integrate with 7 NC points
template<>
template<typename ValueType>
ValueType NC<7>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
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
  
// integrate with 8 NC points
template<>
template<typename ValueType>
ValueType NC<8>::
//computeIntegral(const typename std::array<ValueType, numberEvaluations()>::const_iterator evaluationsIter)
computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluationsIter)
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