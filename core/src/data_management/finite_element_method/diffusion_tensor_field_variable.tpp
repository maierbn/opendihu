#include "data_management/finite_element_method/diffusion_tensor_field_variable.h"

namespace Data
{


template<typename FunctionSpaceType>
DiffusionTensorFieldVariable<FunctionSpaceType>::
DiffusionTensorFieldVariable(PyObject *specificSettings) :
  DiffusionTensorBase<FunctionSpaceType>::DiffusionTensorBase(specificSettings)
{
  LOG(DEBUG) << "construct directional diffusion tensor";
  if (VLOG_IS_ON(1))
  {
    PythonUtility::printDict(this->specificSettings_);
  }

  this->diffusionTensor_ = this->parseDiffusionTensor("diffusionTensor");
}

template<typename FunctionSpaceType>
void DiffusionTensorFieldVariable<FunctionSpaceType>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction, bool useAdditionalDiffusionTensor)
{
  LOG(DEBUG) << "DiffusionTensorFieldVariable::initialize";
  useAdditionalDiffusionTensor_ = useAdditionalDiffusionTensor;
  direction_ = direction;

  if (useAdditionalDiffusionTensor_)
  {
    this->additionalDiffusionTensor_ = this->parseDiffusionTensor("extracellularDiffusionTensor");
  }
}

template<typename FunctionSpaceType>
MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> DiffusionTensorFieldVariable<FunctionSpaceType>::
diffusionTensor(element_no_t elementNoLocal, const std::array<double,FunctionSpaceType::dim()> xi) const
{
  const int D = FunctionSpaceType::dim();

  // get the values at dofs of the element
  std::array<Vec3,FunctionSpaceType::nDofsPerElement()> elementalValues;
  assert(direction_);
  direction_->getElementValues(elementNoLocal, elementalValues);

  // get the function space from field variable
  std::shared_ptr<FunctionSpaceType> functionSpace = direction_->functionSpace();

  // get the interpolated value at xi coordinates inside the element
  Vec3 directionVector = functionSpace->template interpolateValueInElement<3>(elementalValues,xi);

  MathUtility::Matrix<D,D> diffusionTensor = this->diffusionTensor_;

  // if the extracellular diffusion tensor should be added
  if (useAdditionalDiffusionTensor_)
  {
    //LOG(DEBUG) << " add extracellular diffusion tensor " << this->additionalDiffusionTensor_;

    // add extracellular diffusion tensor
    diffusionTensor = diffusionTensor + this->additionalDiffusionTensor_;
  }

  VLOG(3) << "directionVector: " << directionVector;
  //LOG(DEBUG) << "diffusionTensor before rotation: " << diffusionTensor;

  // rotate diffusion tensor in fiber direction
  //MathUtility::rotateMatrix(diffusionTensor, directionVector);

  //LOG(DEBUG) << "diffusionTensor after rotation: " << diffusionTensor;

  return diffusionTensor;
}

}  // namespace
