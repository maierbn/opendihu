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
template<int D, typename GaussQuadrature>
class TriangularPrism :
  public TriangularPrismBase<GaussQuadrature>
{
public:
  //! return the sampling points, i.e. newton-cotes points that are needed for the quadrature. The list may not be in ascending order, but the order matches the order required in integrate
  static std::array<std::array<double,D>,TriangularPrismBase<GaussQuadrature>::numberEvaluations()> samplingPoints();

  //! return the quadrature weights
  static const std::array<double,TriangularPrismBase<GaussQuadrature>::numberEvaluations()> quadratureWeights();

  //! Compute the integral from evaluations at the integration points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same integration points for all.
  //! The number of entries `nEntriesEvaluationArray` of the input array may be larger than the actual number of evaluations, such that
  //! the same array can be used for normal hexahedral elements and triangular prism elements.
  template<typename ValueType, long unsigned int nEntriesEvaluationArray>
  static ValueType computeIntegral(const typename std::array<ValueType,nEntriesEvaluationArray> &evaluations);

  //! Reorder the entries in the evaluations array such that it fits for the triangular prism element, the values were assembled for a 3D hexahedral element
  //! @param usesDofPairs if the evaluations have entries for (i,j)-pairs of dofs, e.g., as used in stiffness and mass matrix. If false, there is only one list of dofs, then isQuadraticElement1 does not matter and nEntriesPerDof1 should be 1
  template<typename ValueType, long unsigned int nEntriesEvaluationArray>
  static void adjustEntriesforPrism(typename std::array<ValueType,nEntriesEvaluationArray> &evaluations,
                                    Mesh::face_or_edge_t edge, bool usesDofPairs, bool isQuadraticElement0, int nEntriesPerDof0,
                                    bool isQuadraticElement1=false, int nEntriesPerDof1=1);
};

/** specialization for D=3, for anything else it does nothing
 */
template<typename GaussQuadrature>
class TriangularPrism<3,GaussQuadrature> :
  public TriangularPrismBase<GaussQuadrature>
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

  //! Reorder the entries in the evaluations array such that it fits for the triangular prism element, the values were assembled for a 3D hexahedral element
  //! @param usesDofPairs if the evaluations have entries for (i,j)-pairs of dofs, e.g., as used in stiffness and mass matrix. If false, there is only one list of dofs, then isQuadraticElement1 does not matter and nEntriesPerDof1 should be 1
  template<typename ValueType, long unsigned int nEntriesEvaluationArray>
  static void adjustEntriesforPrism(typename std::array<ValueType,nEntriesEvaluationArray> &evaluations,
                                    Mesh::face_or_edge_t edge, bool usesDofPairs, bool isQuadraticElement0, int nEntriesPerDof0,
                                    bool isQuadraticElement1=false, int nEntriesPerDof1=1);

};

} // namespace

#include "quadrature/triangular_prism.tpp"
