#include "data_management/finite_element_method/diffusion_tensor_base.h"

namespace Data
{


template<typename FunctionSpaceType>
DiffusionTensorBase<FunctionSpaceType>::
DiffusionTensorBase(PythonConfig specificSettings) :
  specificSettings_(specificSettings)
{
  LOG(DEBUG) << "initialize diffusion tensor";
}

template<typename FunctionSpaceType>
MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> DiffusionTensorBase<FunctionSpaceType>::
parseDiffusionTensor(std::string settingsKey)
{
  const int D = FunctionSpaceType::dim();

  MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> result;
  if (this->specificSettings_.hasKey(settingsKey))
  {
    // create identity matrix as default values
    std::array<double, D*D> defaultValue({0});
    for (int i = 0; i < D; i++)
    {
      defaultValue[i*D] = 1.0;
    }

    // get diffusion tensor from config as array with D*D entries
    result = this->specificSettings_.template getOptionArray<double, D*D>(settingsKey, defaultValue);

    LOG(DEBUG) << "parsed diffusionTensor \"" << settingsKey << "\": " << result;
  }
  else
  {
    LOG(FATAL) << "Diffusion tensor \"" << settingsKey << "\" not found in config.";
  }
  return result;
}


}  // namespace
