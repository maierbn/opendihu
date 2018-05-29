#include "quadrature/clenshaw_curtis.h"

#include <array>
#include <cmath>

namespace Quadrature
{

// 1 ClenshawCurtis point
template<>
std::array<double, 1> ClenshawCurtis<1>::
samplingPoints()
{
  return std::array<double, 1>{0.5};
}

// 2 ClenshawCurtis points
template<>
std::array<double, 2> ClenshawCurtis<2>::
samplingPoints()
{
  return std::array<double, 2>{
    0.,
      1.
  };
}

// 3 ClenshawCurtis points
template<>
std::array<double, 3> ClenshawCurtis<3>::
samplingPoints()
{
  return std::array<double, 3>{
    0.,
      0.5,
      1.
  };
}

// 4 ClenshawCurtis points
template<>
std::array<double, 4> ClenshawCurtis<4>::
samplingPoints()
{
  return std::array<double, 4>{
    0.,
      0.25,
      0.75,
      1.
  };
}

// 5 ClenshawCurtis points
template<>
std::array<double, 5> ClenshawCurtis<5>::
samplingPoints()
{
  return std::array<double, 5>{
    0.,
      (-0.7071067811865475 + 1.) / 2.,
      0.5,
      (0.7071067811865475 + 1.) / 2.,
      1.,
  };
}

// 7 ClenshawCurtis points
template<>
std::array<double, 6> ClenshawCurtis<6>::
samplingPoints()
{
  return std::array<double, 6>{
    0.,
      0.095491502812526274,
      0.34549150281252627,
      0.65450849718747373,
      0.90450849718747373,
      1
  };
}

// 7 ClenshawCurtis points
template<>
std::array<double, 7> ClenshawCurtis<7>::
samplingPoints()
{
  return std::array<double, 7>{
    0.,
      0.066987298107780646,
      0.24999999999999994,
      0.49999999999999994,
      0.74999999999999989,
      0.93301270189221941,
      1
  };
}

// 64 ClenshawCurtis points
template<>
std::array<double, 64> ClenshawCurtis<64>::
samplingPoints()
{
  return std::array<double, 64>{
    0.,
      0.00062153939053882779,
      0.0024846123172992951,
      0.0055845868874357385,
      0.0099137560757280863,
      0.015461356885461019,
      0.022213597106929606,
      0.030153689607045786,
      0.039261894064796188,
      0.049515566048790427,
      0.060889213314885726,
      0.073354559183922141,
      0.086880612842002547,
      0.10143374638853875,
      0.11697777844051099,
      0.13347406408508683,
      0.15088159095696357,
      0.16915708120157025,
      0.1882550990706332,
      0.20812816388260502,
      0.22872686806712028,
      0.24999999999999994,
      0.27189467132341849,
      0.29435644843469416,
      0.31732948781680237,
      0.34075667487415767,
      0.36457976592849739,
      0.38873953302184278,
      0.41317591116653479,
      0.43782814767625738,
      0.46263495320678782,
      0.48753465413096347,
      0.51246534586903636,
      0.53736504679321218,
      0.56217185232374245,
      0.58682408883346515,
      0.61126046697815717,
      0.6354202340715025,
      0.65924332512584216,
      0.68267051218319752,
      0.70564355156530567,
      0.72810532867658151,
      0.74999999999999989,
      0.77127313193287961,
      0.79187183611739487,
      0.81174490092936669,
      0.83084291879842964,
      0.84911840904303637,
      0.866525935914913,
      0.88302222155948895,
      0.89856625361146114,
      0.91311938715799734,
      0.9266454408160778,
      0.93911078668511427,
      0.95048443395120952,
      0.96073810593520381,
      0.96984631039295421,
      0.97778640289307039,
      0.98453864311453887,
      0.99008624392427191,
      0.99441541311256421,
      0.99751538768270076,
      0.99937846060946112,
      1,
  };
}

};  // namespace