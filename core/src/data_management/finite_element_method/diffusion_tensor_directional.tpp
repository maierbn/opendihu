#include "data_management/finite_element_method/diffusion_tensor_directional.h"

namespace Data
{


template<typename FunctionSpaceType>
DiffusionTensorDirectional<FunctionSpaceType>::
DiffusionTensorDirectional(PythonConfig specificSettings) :
  DiffusionTensorBase<FunctionSpaceType>::DiffusionTensorBase(specificSettings)
{
  LOG(DEBUG) << "construct directional diffusion tensor";
  if (VLOG_IS_ON(1))
  {
    PythonUtility::printDict(this->specificSettings_.pyObject());
  }

  this->diffusionTensor_ = this->parseDiffusionTensor("diffusionTensor");
}

template<typename FunctionSpaceType>
void DiffusionTensorDirectional<FunctionSpaceType>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
           bool useAdditionalDiffusionTensor)
{
  LOG(DEBUG) << "DiffusionTensorDirectional::initialize";
  useAdditionalDiffusionTensor_ = useAdditionalDiffusionTensor;
  direction_ = direction;
  spatiallyVaryingPrefactor_ = spatiallyVaryingPrefactor;

  if (useAdditionalDiffusionTensor_)
  {
    this->additionalDiffusionTensor_ = this->parseDiffusionTensor("extracellularDiffusionTensor");
  }
}

template<typename FunctionSpaceType>
MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> DiffusionTensorDirectional<FunctionSpaceType>::
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

  if (!std::isfinite(directionVector[0]) || !std::isfinite(directionVector[1]) || !std::isfinite(directionVector[2]))
    directionVector = Vec3({0.0,0.0,1.0});

  MathUtility::Matrix<D,D> diffusionTensor = this->diffusionTensor_;

  // if the extracellular diffusion tensor should be added
  if (useAdditionalDiffusionTensor_)
  {
    //LOG(DEBUG) << " add extracellular diffusion tensor " << this->additionalDiffusionTensor_;

    // add extracellular diffusion tensor
    diffusionTensor = diffusionTensor + this->additionalDiffusionTensor_;
  }

  //VLOG(3) << "directionVector: " << directionVector;
  //VLOG(2) << "diffusionTensor before rotation: " << diffusionTensor;

  // rotate diffusion tensor in fiber direction
  MathUtility::rotateMatrix(diffusionTensor, directionVector);

  //VLOG(2) << "diffusionTensor after rotation: " << diffusionTensor;

  // if there is a relative factor
  double spatiallyVaryingPrefactor = 1;
  if (spatiallyVaryingPrefactor_)
  {
    std::array<double,FunctionSpaceType::nDofsPerElement()> spatiallyVaryingPrefactorElementalValues;
    spatiallyVaryingPrefactor_->getElementValues(elementNoLocal, spatiallyVaryingPrefactorElementalValues);
    spatiallyVaryingPrefactor = functionSpace->interpolateValueInElement(spatiallyVaryingPrefactorElementalValues, xi);

    diffusionTensor *= spatiallyVaryingPrefactor;
    if (VLOG_IS_ON(1))
      VLOG(1) << "elementNoLocal " << elementNoLocal << ", factor: " << spatiallyVaryingPrefactor << ", scaled diffusionTensor: " << diffusionTensor;
  }

  // check if diffusion tensor contains nan values
  if (D == 3)
  {
    if (!std::isfinite(diffusionTensor[0]) || !std::isfinite(diffusionTensor[1]) || !std::isfinite(diffusionTensor[2])
      || !std::isfinite(diffusionTensor[3]) || !std::isfinite(diffusionTensor[4]) || !std::isfinite(diffusionTensor[5])
      || !std::isfinite(diffusionTensor[6]) || !std::isfinite(diffusionTensor[7]) || !std::isfinite(diffusionTensor[8]))
    {
      LOG(ERROR) << "Directional diffusion tensor contains nans or infs.";
      LOG(INFO) << "elementNoLocal: " << elementNoLocal << ", xi: " << xi;
      LOG(INFO) << "directionVector: " << directionVector << " elemental direction values: " << elementalValues;
      LOG(INFO) << "diffusionTensor: " << diffusionTensor << ", this->additionalDiffusionTensor: " << this->additionalDiffusionTensor_
        << ", spatiallyVaryingPrefactor: " << spatiallyVaryingPrefactor;
     }
  }
  else if (D == 2)
  {
    if (!std::isfinite(diffusionTensor[0]) || !std::isfinite(diffusionTensor[1]) || !std::isfinite(diffusionTensor[2])
      || !std::isfinite(diffusionTensor[3]))
    {
      LOG(ERROR) << "Directional diffusion tensor contains nans or infs.";
      LOG(INFO) << "elementNoLocal: " << elementNoLocal << ", xi: " << xi;
      LOG(INFO) << "directionVector: " << directionVector << " elemental direction values: " << elementalValues;
      LOG(INFO) << "diffusionTensor: " << diffusionTensor << ", this->additionalDiffusionTensor: " << this->additionalDiffusionTensor_
        << ", spatiallyVaryingPrefactor: " << spatiallyVaryingPrefactor;
     }
  }
  else if (D == 1)
  {
    if (!std::isfinite(diffusionTensor[0]))
    {
      LOG(ERROR) << "Directional diffusion tensor contains nans.";
      LOG(INFO) << "elementNoLocal: " << elementNoLocal << ", xi: " << xi;
      LOG(INFO) << "directionVector: " << directionVector << " elemental direction values: " << elementalValues;
      LOG(INFO) << "diffusionTensor: " << diffusionTensor << ", this->additionalDiffusionTensor: " << this->additionalDiffusionTensor_
        << ", spatiallyVaryingPrefactor: " << spatiallyVaryingPrefactor;
     }
  }

  return diffusionTensor;
}

}  // namespace
