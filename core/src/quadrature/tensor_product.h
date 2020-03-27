#pragma once

#include <Python.h>  // has to be the first included header
#include <type_traits>
#include <cmath>
#include <Vc/Vc>

#include "quadrature/quadrature.h"

namespace Quadrature
{

template <unsigned int D, typename Quadrature>
class TensorProductBase
{
public:
  static constexpr int numberEvaluations();
};

/** Integration in D dimension using a tensor product of 1D integrations
 */
template <unsigned int D, typename Quadrature>
class TensorProduct : public TensorProductBase<D, Quadrature>
{
public:
};

// partial specialization for 1D
template<typename Quadrature>
class TensorProduct<1,Quadrature> : public TensorProductBase<1,Quadrature>
{
public:
  //! compute the integral for a scalar value
  template<typename ValueType>
  static ValueType computeIntegral(const Vc::array<ValueType, TensorProductBase<1,Quadrature>::numberEvaluations()> &evaluations);

  //! get the sampling points, i.e. points where the function needs to be evaluated
  static Vc::array<Vc::array<double,1>, TensorProductBase<1,Quadrature>::numberEvaluations()> samplingPoints();
};

// partial specialization for 2D
template<typename Quadrature>
class TensorProduct<2,Quadrature> : public TensorProductBase<2,Quadrature>
{
public:
  //! compute the integral for a scalar value
  template<typename ValueType>
  static ValueType computeIntegral(const Vc::array<ValueType, TensorProductBase<2,Quadrature>::numberEvaluations()> &evaluations);

  //! get the sampling points, i.e. points where the function needs to be evaluated
  static Vc::array<Vc::array<double,2>, TensorProductBase<2,Quadrature>::numberEvaluations()> samplingPoints();
};

// partial specialization for 3D
template<typename Quadrature>
class TensorProduct<3,Quadrature> : public TensorProductBase<3,Quadrature>
{
public:
  //! compute the integral for a scalar value
  template<typename ValueType>
  static ValueType computeIntegral(const Vc::array<ValueType, TensorProductBase<3,Quadrature>::numberEvaluations()> &evaluations);

  //! get the sampling points, i.e. points where the function needs to be evaluated
  static Vc::array<Vc::array<double,3>, TensorProductBase<3,Quadrature>::numberEvaluations()> samplingPoints();
};

} // namespace
#include "quadrature/tensor_product.tpp"
