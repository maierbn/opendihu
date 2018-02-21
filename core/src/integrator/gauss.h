#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <type_traits>

#include "integrator/integrator.h"

namespace Integrator
{

/** Implements Gauss quadrature with NumberGaussPoints. It is capable of exactly integrating polynomials of degree 2*NumberGaussPoints-1.
  */
template<unsigned int NumberGaussPoints>
class Gauss : public Integrator
{
public:
  //! return the number of evaluations that are needed for a 1D quadrature
  static constexpr int numberEvaluations();
  
  //! return the sampling points, i.e. Gauss points that are needed for the quadrature. The list may not be in ascending order, but the order matches the order required in integrate
  static std::array<double, NumberGaussPoints> samplingPoints();
  
  //! Compute the integral from evaluations at the gauss points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same gauss points for all.
  template<typename ValueType>
  static ValueType integrate(const typename std::array<ValueType,numberEvaluations()>::const_iterator evaluations);
  
  //! Compute the integral from evaluations at the gauss points.
  //! If a std::array is given for ValueType, compute separate integrals for each component with the same gauss points for all.
  template<typename ValueType>
  static ValueType integrate(const typename std::array<ValueType,numberEvaluations()> &evaluations);
};
  
};

#include "integrator/gauss.tpp"
