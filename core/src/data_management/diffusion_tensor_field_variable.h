#pragma once

#include <Python.h>  // has to be the first included header

#include "utility/math_utility.h"

namespace Data
{

template<typename FunctionSpaceType>
class DiffusionTensorFieldVariable
{
public:

  //! dummy method that does nothing
  void initialize(){};

  //! initialize diffusion tensor from field variable with directional field
  void initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction, int multidomainNCompartments = 0);

  //! return diffusion tensor
  MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> diffusionTensor(element_no_t elementNoLocal, const std::array<double,FunctionSpaceType::dim()> xi) const;

private:
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction_;   ///< direction of the diffusion tensor
  int multidomainNCompartments_;    ///< if the diffusion tensor should be set as (sigma_i + sigma_e) where sigma_i is multidomainNCompartments_*normal diffusion tensor from directions and sigma_e is in z direction
};

}  // namespace

#include "data_management/diffusion_tensor_field_variable.tpp"
