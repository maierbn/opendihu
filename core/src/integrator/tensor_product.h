#pragma once

#include <Python.h>  // has to be the first included header
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
};

/** Integration in D dimension using a tensor product of 1D integrations
 */
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
  //! compute the integral for a scalar value
  template<typename ValueType>
  static ValueType integrate(const std::array<ValueType, TensorProductBase<1,Integrator>::numberEvaluations()> &evaluations);
  
  //! get the sampling points, i.e. points where the function needs to be evaluated
  static std::array<std::array<double,1>, TensorProductBase<1,Integrator>::numberEvaluations()> samplingPoints();
};

// partial specialization for 2D
template<typename Integrator>
class TensorProduct<2,Integrator> : public TensorProductBase<2,Integrator>
{
public:
  //! compute the integral for a scalar value
  template<typename ValueType>
  static ValueType integrate(const std::array<ValueType, TensorProductBase<2,Integrator>::numberEvaluations()> &evaluations);
  
  //! get the sampling points, i.e. points where the function needs to be evaluated
  static std::array<std::array<double,2>, TensorProductBase<2,Integrator>::numberEvaluations()> samplingPoints();
};

// partial specialization for 3D
template<typename Integrator>
class TensorProduct<3,Integrator> : public TensorProductBase<3,Integrator>
{
public:
  //! compute the integral for a scalar value
  template<typename ValueType>
  static ValueType integrate(const std::array<ValueType, TensorProductBase<3,Integrator>::numberEvaluations()> &evaluations);
  
  //! get the sampling points, i.e. points where the function needs to be evaluated
  static std::array<std::array<double,3>, TensorProductBase<3,Integrator>::numberEvaluations()> samplingPoints();
};

};  // namespace
#include "integrator/tensor_product.tpp"
