#include "quadrature/triangular_prism.h"

#include <array>
#include <cmath>

#include "quadrature/gauss.h"

namespace Quadrature
{

// --------------------
// sampling points
// 1 point in the triangle, 1 point in z direction
template<>
std::array<Vec3,1> TriangularPrism<3,Gauss<1>>::
samplingPoints()
{
  return std::array<Vec3,1>{
    Vec3{1/3., 1/3., Gauss<1>::samplingPoints()[0]}
  };
}

// 3 points in the triangle, 2 points in z direction
template<>
std::array<Vec3,6> TriangularPrism<3,Gauss<2>>::
samplingPoints()
{
  return std::array<Vec3,6>{
    Vec3{0.5, 0,   Gauss<2>::samplingPoints()[0]},
    Vec3{0.5, 0.5, Gauss<2>::samplingPoints()[0]},
    Vec3{0,   0.5, Gauss<2>::samplingPoints()[0]},
    Vec3{0.5, 0,   Gauss<2>::samplingPoints()[1]},
    Vec3{0.5, 0.5, Gauss<2>::samplingPoints()[1]},
    Vec3{0,   0.5, Gauss<2>::samplingPoints()[1]}
  };
}

// 3 points in the triangle, 2 points in z direction
template<>
std::array<Vec3,12> TriangularPrism<3,Gauss<3>>::
samplingPoints()
{
  return std::array<Vec3,12>{
    Vec3{1/3., 1/3., Gauss<3>::samplingPoints()[0]},
    Vec3{0.6,  0.2,  Gauss<3>::samplingPoints()[0]},
    Vec3{0.2,  0.2,  Gauss<3>::samplingPoints()[0]},
    Vec3{0.2,  0.6,  Gauss<3>::samplingPoints()[0]},
    Vec3{1/3., 1/3., Gauss<3>::samplingPoints()[1]},
    Vec3{0.6,  0.2,  Gauss<3>::samplingPoints()[1]},
    Vec3{0.2,  0.2,  Gauss<3>::samplingPoints()[1]},
    Vec3{0.2,  0.6,  Gauss<3>::samplingPoints()[1]},
    Vec3{1/3., 1/3., Gauss<3>::samplingPoints()[2]},
    Vec3{0.6,  0.2,  Gauss<3>::samplingPoints()[2]},
    Vec3{0.2,  0.2,  Gauss<3>::samplingPoints()[2]},
    Vec3{0.2,  0.6,  Gauss<3>::samplingPoints()[2]}
  };
}

// --------------------
// quadrature weights

// 1 point in the triangle, 1 point in z direction
template<>
const std::array<double,1> TriangularPrism<3,Gauss<1>>::
quadratureWeights()
{
  return std::array<double,1>{
    0.5 * 1.0 * Gauss<1>::quadratureWeights()[0]
  };
}

// 3 points in the triangle, 2 points in z direction
template<>
const std::array<double,6> TriangularPrism<3,Gauss<2>>::
quadratureWeights()
{
  return std::array<double,6>{
    0.5 * 1/3. * Gauss<2>::quadratureWeights()[0],
    0.5 * 1/3. * Gauss<2>::quadratureWeights()[0],
    0.5 * 1/3. * Gauss<2>::quadratureWeights()[0],
    0.5 * 1/3. * Gauss<2>::quadratureWeights()[1],
    0.5 * 1/3. * Gauss<2>::quadratureWeights()[1],
    0.5 * 1/3. * Gauss<2>::quadratureWeights()[1]
  };
}

// 3 points in the triangle, 2 points in z direction
template<>
const std::array<double,12> TriangularPrism<3,Gauss<3>>::
quadratureWeights()
{
  return std::array<double,12>{
    0.5 * -27/48. * Gauss<3>::quadratureWeights()[0],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[0],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[0],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[0],
    0.5 * -27/48. * Gauss<3>::quadratureWeights()[1],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[1],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[1],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[1],
    0.5 * -27/48. * Gauss<3>::quadratureWeights()[2],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[2],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[2],
    0.5 * 25/48. * Gauss<3>::quadratureWeights()[2]
  };
}

}  // namespace
