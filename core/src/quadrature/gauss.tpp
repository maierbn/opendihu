#include "quadrature/gauss.h"

#include <array>

namespace Quadrature
{

template<unsigned int NumberGaussPoints>
constexpr int Gauss<NumberGaussPoints>::numberEvaluations()
{
  return NumberGaussPoints;
}

template<unsigned int NumberGaussPoints>
template<typename ValueType>
ValueType Gauss<NumberGaussPoints>::computeIntegral(const typename std::array<ValueType,Gauss<NumberGaussPoints>::numberEvaluations()> &evaluations)
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

// quadrature with 2 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<2>::
computeIntegral(const typename std::array<ValueType,2>::const_iterator evaluationsIter)
{
  ValueType result{};
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
computeIntegral(const typename std::array<ValueType,3>::const_iterator evaluationsIter)
{
  ValueType result{};
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
computeIntegral(const typename std::array<ValueType,4>::const_iterator evaluationsIter)
{
  ValueType result{};
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
computeIntegral(const typename std::array<ValueType,5>::const_iterator evaluationsIter)
{
  ValueType result{};
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
computeIntegral(const typename std::array<ValueType,6>::const_iterator evaluationsIter)
{
  ValueType result{};
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
computeIntegral(const typename std::array<ValueType,7>::const_iterator evaluationsIter)
{
  ValueType result{};
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

// integrate with 8 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<8>::
computeIntegral(const typename std::array<ValueType,8>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,8> weights {
    0.3626837833783620/2.,
    0.3626837833783620/2.,  
    0.3137066458778873/2.,  
    0.3137066458778873/2.,  
    0.2223810344533745/2.,  
    0.2223810344533745/2.,  
    0.1012285362903763/2.,  
    0.1012285362903763/2.
  };
  for (int i = 0; i < 8; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}


// integrate with 10 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<10>::
computeIntegral(const typename std::array<ValueType,10>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,10> weights {
    0.2955242247147529/2.,  
    0.2955242247147529/2.,  
    0.2692667193099963/2.,  
    0.2692667193099963/2.,  
    0.2190863625159820/2.,  
    0.2190863625159820/2.,  
    0.1494513491505806/2.,  
    0.1494513491505806/2.,  
    0.0666713443086881/2.,  
    0.0666713443086881/2.
  };
  for (int i = 0; i < 10; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}


// integrate with 12 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<12>::
computeIntegral(const typename std::array<ValueType,12>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,12> weights {
    0.2491470458134028/2.,  
    0.2491470458134028/2.,  
    0.2334925365383548/2.,  
    0.2334925365383548/2.,  
    0.2031674267230659/2.,  
    0.2031674267230659/2.,  
    0.1600783285433462/2.,  
    0.1600783285433462/2.,  
    0.1069393259953184/2.,  
    0.1069393259953184/2.,
    0.0471753363865118/2.,
    0.0471753363865118/2.
  };
  for (int i = 0; i < 12; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}


// integrate with 16 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<16>::
computeIntegral(const typename std::array<ValueType,16>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,16> weights {
    0.1894506104550685/2.,  
    0.1894506104550685/2.,  
    0.1826034150449236/2.,  
    0.1826034150449236/2.,  
    0.1691565193950025/2.,  
    0.1691565193950025/2.,  
    0.1495959888165767/2.,  
    0.1495959888165767/2.,  
    0.1246289712555339/2.,  
    0.1246289712555339/2.,
    0.0951585116824928/2.,
    0.0951585116824928/2.,
    0.0622535239386479/2.,
    0.0622535239386479/2.,
    0.0271524594117541/2.,
    0.0271524594117541/2.
  };
  for (int i = 0; i < 16; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}


// integrate with 20 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<20>::
computeIntegral(const typename std::array<ValueType,20>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,20> weights {
    0.1527533871307258/2.,  
    0.1527533871307258/2.,  
    0.1491729864726037/2.,  
    0.1491729864726037/2.,  
    0.1420961093183820/2.,  
    0.1420961093183820/2.,  
    0.1316886384491766/2.,  
    0.1316886384491766/2.,  
    0.1181945319615184/2.,  
    0.1181945319615184/2.,
    0.1019301198172404/2.,
    0.1019301198172404/2.,
    0.0832767415767048/2.,
    0.0832767415767048/2.,
    0.0626720483341091/2.,
    0.0626720483341091/2.,
    0.0406014298003869/2.,
    0.0406014298003869/2.,
    0.0176140071391521/2.,
    0.0176140071391521/2.
  };
  for (int i = 0; i < 20; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}


// integrate with 24 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<24>::
computeIntegral(const typename std::array<ValueType,24>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,24> weights {
    0.1279381953467522/2.,  
    0.1279381953467522/2.,
    0.1258374563468283/2.,  
    0.1258374563468283/2.,  
    0.1216704729278034/2.,  
    0.1216704729278034/2.,  
    0.1155056680537256/2.,  
    0.1155056680537256/2.,  
    0.1074442701159656/2.,  
    0.1074442701159656/2.,
    0.0976186521041139/2.,
    0.0976186521041139/2.,
    0.0861901615319533/2.,
    0.0861901615319533/2.,
    0.0733464814110803/2.,
    0.0733464814110803/2.,
    0.0592985849154368/2.,
    0.0592985849154368/2.,
    0.0442774388174198/2.,
    0.0442774388174198/2.,
    0.0285313886289337/2.,
    0.0285313886289337/2.,
    0.0123412297999872/2.,
    0.0123412297999872/2.
  };
  for (int i = 0; i < 24; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

// integrate with 64 Gauss points
template<>
template<typename ValueType>
ValueType Gauss<64>::
computeIntegral(const typename std::array<ValueType,64>::const_iterator evaluationsIter)
{
  ValueType result{};
  std::array<double,64> weights {
    0.0486909570091397/2.,
    0.0486909570091397/2.,
    0.0485754674415034/2.,
    0.0485754674415034/2.,
    0.0483447622348030/2.,
    0.0483447622348030/2.,
    0.0479993885964583/2.,
    0.0479993885964583/2.,
    0.0475401657148303/2.,
    0.0475401657148303/2.,
    0.0469681828162100/2.,
    0.0469681828162100/2.,
    0.0462847965813144/2.,
    0.0462847965813144/2.,
    0.0454916279274181/2.,
    0.0454916279274181/2.,
    0.0445905581637566/2.,
    0.0445905581637566/2.,
    0.0435837245293235/2.,
    0.0435837245293235/2.,
    0.0424735151236536/2.,
    0.0424735151236536/2.,
    0.0412625632426235/2.,
    0.0412625632426235/2.,
    0.0399537411327203/2.,
    0.0399537411327203/2.,
    0.0385501531786156/2.,
    0.0385501531786156/2.,
    0.0370551285402400/2.,
    0.0370551285402400/2.,
    0.0354722132568824/2.,
    0.0354722132568824/2.,
    0.0338051618371416/2.,
    0.0338051618371416/2.,
    0.0320579283548516/2.,
    0.0320579283548516/2.,
    0.0302346570724025/2.,
    0.0302346570724025/2.,
    0.0283396726142595/2.,
    0.0283396726142595/2.,
    0.0263774697150547/2.,
    0.0263774697150547/2.,
    0.0243527025687109/2.,
    0.0243527025687109/2.,
    0.0222701738083833/2.,
    0.0222701738083833/2.,
    0.0201348231535302/2.,
    0.0201348231535302/2.,
    0.0179517157756973/2.,
    0.0179517157756973/2.,
    0.0157260304760247/2.,
    0.0157260304760247/2.,
    0.0134630478967186/2.,
    0.0134630478967186/2.,
    0.0111681394601311/2.,
    0.0111681394601311/2.,
    0.0088467598263639/2.,
    0.0088467598263639/2.,
    0.0065044579689784/2.,
    0.0065044579689784/2.,
    0.0041470332605625/2.,
    0.0041470332605625/2.,
    0.0017832807216964/2.,
    0.0017832807216964/2.
  };
  for (int i = 0; i < 64; i++)
  {
    result += weights[i] * (*(evaluationsIter+i));
  }
  return result;
}

}  // namespace
