#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <type_traits>

#include "quadrature/quadrature.h"
#include "control/types.h"
#include "mesh/face_or_edge_t.h"
#include "quadrature/triangular_prism_base.h"

namespace Quadrature
{

/** Quadrature of a triangle in x-y plane and a rectangle in z direction (i.e., a triangular prism).
 *  The sampling points of the triangle are taken from Zienkiewicz, Taylor (2005) "The Finite Element Method: Its Basis and Fundamentals", Table 5.3.
*/
template<typename GaussQuadrature>
class TriangularPrism : public TriangularPrismBase<GaussQuadrature>
{
public:
  //! return the sampling points, i.e. newton-cotes points that are needed for the quadrature. The list may not be in ascending order, but the order matches the order required in integrate
  static std::array<Vec3,TriangularPrismBase<GaussQuadrature>::numberEvaluations()> samplingPoints();

  //! return the quadrature weights
  static const std::array<double,TriangularPrismBase<GaussQuadrature>::numberEvaluations()> quadratureWeights();

  //! Compute the integral from evaluations at the integration points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same integration points for all.
  //! The number of entries `nEntriesEvaluationArray` of the input array may be larger than the actual number of evaluations, such that
  //! the same array can be used for normal hexahedral elements and triangular prism elements.
  template<typename ValueType, long unsigned int nEntriesEvaluationArray>
  static ValueType computeIntegral(const typename std::array<ValueType,nEntriesEvaluationArray> &evaluations);

  //! Reorder the entries in the evaluation array that were assembled for a 3D hexahedral element, such that it fits for the triangular hex element
  template<typename ValueType, long unsigned int nEntriesEvaluationArray>
  static void adjustEntriesforPrism(typename std::array<ValueType,nEntriesEvaluationArray> &evaluations,
                                    Mesh::face_or_edge_t edge, bool isQuadraticElement0, int nEntriesPerDof0,
                                    bool isQuadraticElement1, int nEntriesPerDof1);

};

} // namespace

#include "quadrature/triangular_prism.tpp"
