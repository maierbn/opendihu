#include "data_management/diffusion_tensor.h"

namespace Data
{


template <int D>
void DiffusionTensor<D>::
initialize(PyObject *settings)
{
  LOG(DEBUG) << "initialize Diffusion tensor";
  PythonUtility::printDict(settings);
  if (PythonUtility::hasKey(settings, "diffusionTensor"))
  {
    // create identity matrix as default values
    std::array<double, D*D> defaultValue({0});
    for (int i=0; i<D; i++)
      defaultValue[i*D] = 1.0;

    // get diffusion tensor from config as array with D*D entries
    this->diffusionTensor_ = PythonUtility::template getOptionArray<double, D*D>(settings, "diffusionTensor", defaultValue);

    LOG(DEBUG) << "diffusionTensor found: " << this->diffusionTensor_;
  }
  else
  {
    LOG(DEBUG) << "diffusionTensor not found";
    /*
    // create identity matrix as default values
    this->diffusionTensor_ = std::array<double, D*D>({0});
    for (int i=0; i<D; i++)
      this->diffusionTensor_[i*D] = 1.0;
    */
  }
}

template <int D>
const MathUtility::Matrix<D,D> &DiffusionTensor<D>::
diffusionTensor() const
{
  return this->diffusionTensor_;
}

}  // namespace