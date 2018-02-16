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
  static constexpr int numberEvaluations();
  static std::array<double, NumberGaussPoints> samplingPoints();
  static double integrate(std::array<double, numberEvaluations()> &evaluations);
  static double integrate(double *evaluations);

};
  
};

#include "integrator/gauss.tpp"
