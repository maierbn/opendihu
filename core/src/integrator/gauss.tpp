#include "integrator/gauss.h"

#include <array>

namespace Integrator
{

template<unsigned int NumberGaussPoints>
constexpr int Gauss<NumberGaussPoints>::numberEvaluations()
{
  return NumberGaussPoints;
}

template<unsigned int NumberGaussPoints>
double Gauss<NumberGaussPoints>::integrate(std::array<double, numberEvaluations()> &evaluations)
{
  return Gauss<NumberGaussPoints>::integrate(evaluations.data());
}

template<unsigned int NumberGaussPoints>
double Gauss<NumberGaussPoints>::integrate(double *evaluations)
{
  double result = 0.0;
  return result;
}

};