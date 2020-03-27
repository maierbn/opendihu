#pragma once

#include <Python.h>  // has to be the first included header
#include <Vc/Vc>
#include <type_traits>

#include "quadrature/quadrature.h"

namespace Quadrature
{

/** Implements Newton Cotes quadrature with NumberIntegrationPoints. 
 * It is capable of exactly integrating polynomials of degree int((NumberIntegrationPoints+1)/2)*2-1
 * @author: Dominik Sellenthin, SGS Seminar, 2018 
 *              integration points | exact integration of polynomials
 *   rectangle: 1                    1
 *   trapezoid: 2                    1
 *   Simpson:   3                    3
 *   3/8-rule:  4                    3
 *              5                    5 etc.
*/
template<unsigned int NumberIntegrationPoints>
class NewtonCotes : public Quadrature
{
public:
  typedef NewtonCotes<NumberIntegrationPoints> HighOrderQuadrature;   ///< this defines the own class, to be able to generalize code to mixed quadrature

  //! return the number of evaluations that are needed for a 1D quadrature
  static constexpr int numberEvaluations();

  //! return the sampling points, i.e. newton-cotes points that are needed for the quadrature. The list may not be in ascending order, but the order matches the order required in integrate
  static Vc::array<double, NumberIntegrationPoints> samplingPoints();

  //! return the quadrature weights
  static const Vc::array<double, NumberIntegrationPoints> quadratureWeights();

  //! Compute the integral from evaluations at the integration points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same integration points for all.
  template<typename ValueType>
  static ValueType computeIntegral(const typename Vc::array<ValueType,NewtonCotes<NumberIntegrationPoints>::numberEvaluations()>::const_iterator evaluations);
  
  //! Compute the integral from evaluations at the integration points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same integration points for all.
  template<typename ValueType>
  static ValueType computeIntegral(const typename Vc::array<ValueType,NewtonCotes<NumberIntegrationPoints>::numberEvaluations()> &evaluations);
};

} // namespace

#include "quadrature/newton_cotes.tpp"
