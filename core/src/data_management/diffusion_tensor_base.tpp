#include "data_management/diffusion_tensor_base.h"

namespace Data
{


template<typename FunctionSpaceType>
DiffusionTensorBase<FunctionSpaceType>::
DiffusionTensorBase(PyObject *specificSettings) :
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
  if (PythonUtility::hasKey(this->specificSettings_, settingsKey))
  {
    // create identity matrix as default values
    std::array<double, D*D> defaultValue({0});
    for (int i = 0; i < D; i++)
    {
      defaultValue[i*D] = 1.0;
    }

    // get diffusion tensor from config as array with D*D entries
    result = PythonUtility::template getOptionArray<double, D*D>(this->specificSettings_, settingsKey, defaultValue);

    LOG(DEBUG) << "parsed diffusionTensor \"" << settingsKey << "\": " << result;
  }
  else
  {
    LOG(FATAL) << "Diffusion tensor \"" << settingsKey << "\" not found in config.";
  }
  return result;
}


}  // namespace
