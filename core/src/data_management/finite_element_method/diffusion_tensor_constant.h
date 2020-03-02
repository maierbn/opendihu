#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/finite_element_method/diffusion_tensor_base.h"
#include "utility/math_utility.h"
#include "control/python_config/spatial_parameter.h"

namespace Data
{

template<typename FunctionSpaceType>
class DiffusionTensorConstant :
  public DiffusionTensorBase<FunctionSpaceType>
{
public:

  //! constructor
  DiffusionTensorConstant(PythonConfig specificSettings);

  //! read values of diffusion tensor from config
  void initialize(std::shared_ptr<FunctionSpaceType> functionSpace);

  //! return diffusion tensor
  const MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> &diffusionTensor(element_no_t elementNoLocal, const std::array<double,FunctionSpaceType::dim()> xi) const;

private:
  SpatialParameter<FunctionSpaceType,MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()>> diffusionTensor_;  ///< the diffusion/conductivity tensor A in an equation ∇•A∇ = f
};

}  // namespace

#include "data_management/finite_element_method/diffusion_tensor_constant.tpp"
