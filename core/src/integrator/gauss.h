#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <type_traits>

#include "integrator/integrator.h"

namespace Integrator
{

template<unsigned int NumberGaussPoints>
class Gauss : Integrator
{
public:
  static constexpr int numberEvaluations();
  static std::array<double, NumberGaussPoints> samplingPoints();
  static double integrate(std::array<double, numberEvaluations()> &evaluations);
  static double integrate(double *evaluations);

};
  
};

#include "integrator/gauss.tpp"
