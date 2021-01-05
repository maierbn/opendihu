#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <type_traits>

#include "quadrature/quadrature.h"
#include "control/types.h"

namespace Quadrature
{

/** Quadrature of a triangle in x-y plane and a rectangle in z direction (i.e., a triangular prism).
 *  The sampling points of the triangle are taken from Zienkiewicz, Taylor (2005) "The Finite Element Method: Its Basis and Fundamentals", Table 5.3.
*/
template<typename GaussQuadrature>
class TriangularPrismBase : public Quadrature
{
public:
  //! return the number of evaluations that are needed for the quadrature
  static constexpr int numberEvaluations()
  {
    constexpr int n = GaussQuadrature::numberEvaluations();
    constexpr int nTriangle = int(n/2) + n;     // 1 -> 0,  2 -> 1,  3 -> 1
    return nTriangle * n;
  }
};

} // namespace
