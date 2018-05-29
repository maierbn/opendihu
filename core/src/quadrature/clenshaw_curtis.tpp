// weights are calculated corresponding to the implementation in ClenshawCurtisweights()
// which is located in examples/quadrature/dominik/src/quadrature.cpp

#include "quadrature/clenshaw_curtis.h"

#include <array>

namespace Quadrature
{

template<unsigned int NumberIntegrationPoints>
constexpr int ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()
{
  return NumberIntegrationPoints;
}

template<unsigned int NumberIntegrationPoints>
template<typename ValueType>
ValueType ClenshawCurtis<NumberIntegrationPoints>::computeIntegral(const typename std::array<ValueType,ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()> &evaluations)
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
  result = 2. * (*evaluationsIter) * 0.5;
  return result;
}
  
// quadrature with 2 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<2>::
computeIntegral(const typename std::array<ValueType,2>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,2> weights {1., 1.};
  for (int i = 0; i < 2; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}

// quadrature with 3 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<3>::
computeIntegral(const typename std::array<ValueType,3>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,3> weights {1./3., 4./3., 1./3.};
  for (int i = 0; i < 3; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}

// integrate with 4 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<4>::
computeIntegral(const typename std::array<ValueType,4>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,4> weights {
    1./9., 8./9., 8./9., 1./9.
  };
  for (int i = 0; i < 4; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}

// integrate with 5 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<5>::
computeIntegral(const typename std::array<ValueType,5>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,5> weights {
    0.0666666666666666, 0.5333333333333333, 0.7999999999999999, 0.5333333333333333, 0.0666666666666666
  };
  for (int i = 0; i < 5; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}

// integrate with 6 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<6>::
computeIntegral(const typename std::array<ValueType,6>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,6> weights {
    0.040000000000000008,
	0.36074304120001122,
	0.59925695879998886,
	0.59925695879998875,
	0.36074304120001133,
	0.040000000000000008
  };
  for (int i = 0; i < 6; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}

// integrate with 7 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<7>::
computeIntegral(const typename std::array<ValueType,7>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,7> weights {
    0.028571428571428577,
    0.2539682539682539,
    0.45714285714285713,
    0.52063492063492056,
    0.4571428571428573,
    0.2539682539682539,
    0.028571428571428577
  };
  for (int i = 0; i < 7; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}
  
// integrate with 64 ClenshawCurtis points
template<>
template<typename ValueType>
ValueType ClenshawCurtis<64>::
computeIntegral(const typename std::array<ValueType,64>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,64> weights {
    0.00025195263290501448,
    0.0024268255770913318,
    0.0049856828009972485,
    0.0074221088952156119,
    0.0098867917339801568,
    0.012301075649320037,
    0.01470112902792591,
    0.017053372663619865,
    0.01937142575637869,
    0.021635087608271566,
    0.02384983935834243,
    0.026001385652925597,
    0.028091483494470625,
    0.030109091643713892,
    0.032054068367475069,
    0.033917467729185892,
    0.035698154368101788,
    0.037388704313403007,
    0.038987496746689054,
    0.040488323162672368,
    0.041889389342280314,
    0.04318553037633184,	
    0.044374982977791189,
    0.045453526843574302,
    0.046419569176223564,
    0.047269776669786039,
    0.048002824246101541,
    0.048616232217200656,
    0.049109010480342288,
    0.049479514061046455,
    0.049727132112596689,
    0.049851044314039916,
    0.049851044314039916,
    0.049727132112596689,
    0.049479514061046455,
    0.049109010480342288,
    0.048616232217200663,
    0.048002824246101548,
    0.047269776669786052,
    0.046419569176223564,
    0.045453526843574309,
    0.044374982977791223,
    0.04318553037633184,
    0.041889389342280314,
    0.040488323162672389,
    0.038987496746689033,
    0.037388704313403021,
    0.035698154368101795,
    0.033917467729185878,
    0.032054068367475076,
    0.030109091643713896,
    0.028091483494470638,
    0.026001385652925594,
    0.02384983935834243,
    0.021635087608271566,
    0.019371425756378676,
    0.017053372663619868,
    0.014701129027925933,
    0.012301075649320059,
    0.0098867917339801534,
    0.0074221088952156266,
    0.0049856828009972511,
    0.0024268255770913526,
    0.00025195263290501448,
  };
  for (int i = 0; i < 64; i++)
  {
    result += weights[i] * (*(evaluationsIter+i)) * 0.5;
  }
  return result;
}
  
};