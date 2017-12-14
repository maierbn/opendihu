#pragma once

#include <type_traits>
#include <cmath>
#include <array>

#include "integrator/integrator.h"

namespace Integrator 
{

template <unsigned int D, typename Integrator> 
class TensorProductBase
{
public:
  static constexpr int numberEvaluations();
  static constexpr int samplingArraySize();
};

template <unsigned int D, typename Integrator>
class TensorProduct : public TensorProductBase<D, Integrator>
{
public:  
};
  
// partial specialization for 1D
template<typename Integrator>
class TensorProduct<1,Integrator> : public TensorProductBase<1,Integrator>
{
public:
  static double integrate(std::array<double, TensorProductBase<1,Integrator>::numberEvaluations()> &evaluations);
  static std::array<double, TensorProductBase<1,Integrator>::samplingArraySize()> samplingPoints();
};

// partial specialization for 2D
template<typename Integrator>
class TensorProduct<2,Integrator> : public TensorProductBase<2,Integrator>
{
public:
  static double integrate(std::array<double, TensorProductBase<2,Integrator>::numberEvaluations()> &evaluations);
  static std::array<double, TensorProductBase<2,Integrator>::samplingArraySize()> samplingPoints();
};

// partial specialization for 3D
template<typename Integrator>
class TensorProduct<3,Integrator> : public TensorProductBase<3,Integrator>
{
public:
  static double integrate(std::array<double, TensorProductBase<3,Integrator>::numberEvaluations()> &evaluations);
  static std::array<double, TensorProductBase<3,Integrator>::samplingArraySize()> samplingPoints();
};

};  // namespace
#include "integrator/tensor_product.tpp"