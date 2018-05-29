#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <type_traits>

#include "quadrature/integrator.h"

namespace Quadrature
{

/** Implements Clenshaw Curtis quadrature with NumberIntegrationPoints. 
 * It is capable of exactly integrating polynomials of degree NumberIntegrationPoints-1.
 * @author: Dominik Sellenthin, SGS Seminar, 2018
*/
template<unsigned int NumberIntegrationPoints>
class ClenshawCurtis : public Quadrature
{
public:
  typedef ClenshawCurtis<NumberIntegrationPoints> HighOrderQuadrature;   ///< this defines the own class, to be able to generalize code to mixed quadrature

  //! return the number of evaluations that are needed for a 1D quadrature
  static constexpr int numberEvaluations();

  //! return the sampling points, i.e. clenshaw-curtis points that are needed for the quadrature. The list may not be in ascending order, but the order matches the order required in integrate
  static std::array<double, NumberIntegrationPoints> samplingPoints();

  //! Compute the integral from evaluations at the integration points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same integration points for all.
  template<typename ValueType>
  static ValueType computeIntegral(const typename std::array<ValueType,ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()>::const_iterator evaluations);
  
  //! Compute the integral from evaluations at the integration points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same integration points for all.
  template<typename ValueType>
  static ValueType computeIntegral(const typename std::array<ValueType,ClenshawCurtis<NumberIntegrationPoints>::numberEvaluations()> &evaluations);
};

};   // namespace

#include "quadrature/clenshaw_curtis.tpp"