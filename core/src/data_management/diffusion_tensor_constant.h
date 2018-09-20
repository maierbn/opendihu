#pragma once

#include <Python.h>  // has to be the first included header

#include "utility/math_utility.h"

namespace Data
{

template <int D>
class DiffusionTensorConstant
{
public:

  //! read values of diffusion tensor from config
  void initialize(PyObject *settings);

  //! return diffusion tensor
  const MathUtility::Matrix<D,D> &diffusionTensor(element_no_t elementNoLocal, const std::array<double,D> xi) const;

private:
  MathUtility::Matrix<D,D> diffusionTensor_;  ///< the diffusion/conductivity tensor A in an equation ∇•A∇ = f
};

}  // namespace

#include "data_management/diffusion_tensor_constant.tpp"
