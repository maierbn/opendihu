#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/finite_element_method/diffusion_tensor_base.h"
#include "utility/math_utility.h"
#include "control/python_config/spatial_parameter.h"

namespace Data
{

/** A Diffusion tensor that points in a given fiber direction and can have a spatially varying prefactor, i.e.
 *  c(x)*A(f(x)),  f(x)...fibre direction at x, c(x)...scalar prefactor
 *
 **/
template<typename FunctionSpaceType>
class DiffusionTensorDirectional :
  public DiffusionTensorBase<FunctionSpaceType>
{
public:

  //! constructor
  DiffusionTensorDirectional(PythonConfig specificSettings);

  //! initialize diffusion tensor from field variable with directional field, spatiallyVaryingPrefactor can be nullptr, then the factor is constant 1
  void initialize(std::shared_ptr<FunctionSpaceType> functionSpace,
                  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
                  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
                  bool useAdditionalDiffusionTensor);

  //! return diffusion tensor
  MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> diffusionTensor(element_no_t elementNoLocal, const std::array<double,FunctionSpaceType::dim()> xi) const;

private:
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction_;                                  //< direction of the diffusion tensor
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor_;                  //< spatially varying factor with which diffusion tensor will be multiplied
  int multidomainNCompartments_;        //< if the diffusion tensor should be set as (sigma_i + sigma_e) where sigma_i is multidomainNCompartments_*normal diffusion tensor from directions and sigma_e is in z direction

  SpatialParameter<FunctionSpaceType,MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()>> diffusionTensor_;      //< the diffusion/conductivity tensor, such that (1,0,0) is the fiber direction
  SpatialParameter<FunctionSpaceType,MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()>> additionalDiffusionTensor_;   //< an additional diffusion tensor, in fiber direction
  bool useAdditionalDiffusionTensor_;   //< if the additionalDiffusionTensor_ should be added to diffusion tensor
};

}  // namespace

#include "data_management/finite_element_method/diffusion_tensor_directional.tpp"
