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
}

template<typename FunctionSpaceType>
void DiffusionTensorDirectional<FunctionSpaceType>::
initialize(std::shared_ptr<FunctionSpaceType> functionSpace,
           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
           bool useAdditionalDiffusionTensor)
{
  LOG(DEBUG) << "DiffusionTensorDirectional::initialize";
  useAdditionalDiffusionTensor_ = useAdditionalDiffusionTensor;
  direction_ = direction;
  spatiallyVaryingPrefactor_ = spatiallyVaryingPrefactor;
  this->dataFunctionSpace_ = functionSpace;

  // initialize diffusion tensor
  const int D = FunctionSpaceType::dim();

  // create identity matrix as default values
  MathUtility::Matrix<D,D> defaultValue({0});

  for (int i = 0; i < D; i++)
    defaultValue(i,i) = 1.0;

  // parse diffusion tensor
  this->diffusionTensor_.initialize(this->specificSettings_, "diffusionTensor", defaultValue, this->dataFunctionSpace_);

  // parse additional diffusion tensor
  if (useAdditionalDiffusionTensor_)
  {
    const int D = FunctionSpaceType::dim();

    // create identity matrix as default values
    MathUtility::Matrix<D,D> defaultValue({0});

    for (int i = 0; i < D; i++)
      defaultValue(i,i) = 1.0;

    this->additionalDiffusionTensor_.initialize(this->specificSettings_, "extracellularDiffusionTensor", defaultValue, this->dataFunctionSpace_);
  }
}

template<typename FunctionSpaceType>
template<typename double_v_t, typename element_no_v_t>
MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim(),double_v_t> DiffusionTensorDirectional<FunctionSpaceType>::
diffusionTensor(element_no_v_t elementNoLocal, const std::array<double,FunctionSpaceType::dim()> xi) const
{
  const int D = FunctionSpaceType::dim();

  // get the values at dofs of the element
  std::array<VecD<3,double_v_t>,FunctionSpaceType::nDofsPerElement()> elementalValues;
  assert(direction_);
  direction_->getElementValues(elementNoLocal, elementalValues);

  // get the function space from field variable
  std::shared_ptr<FunctionSpaceType> functionSpace = direction_->functionSpace();

  // get the interpolated value at xi coordinates inside the element
  VecD<3,double_v_t> directionVector = functionSpace->template interpolateValueInElement<3>(elementalValues, xi);

  if (!MathUtility::isFinite(directionVector[0]) || !MathUtility::isFinite(directionVector[1]) || !MathUtility::isFinite(directionVector[2]))
  {
    directionVector = VecD<3,double_v_t>({0.0,0.0,1.0});
  }

  MathUtility::Matrix<D,D,double_v_t> diffusionTensor = this->diffusionTensor_.value(elementNoLocal);

  // if the extracellular diffusion tensor should be added
  if (useAdditionalDiffusionTensor_)
  {
    //LOG(DEBUG) << " add extracellular diffusion tensor " << this->additionalDiffusionTensor_;

    // add extracellular diffusion tensor
    diffusionTensor = diffusionTensor + this->additionalDiffusionTensor_.value(elementNoLocal);
  }

  //VLOG(3) << "directionVector: " << directionVector;
  //VLOG(2) << "diffusionTensor before rotation: " << diffusionTensor;

  // rotate diffusion tensor in fiber direction
  MathUtility::rotateMatrix(diffusionTensor, directionVector);

  //VLOG(2) << "diffusionTensor after rotation: " << diffusionTensor;

  // if there is a relative factor
  double_v_t spatiallyVaryingPrefactor = 1;
  if (spatiallyVaryingPrefactor_)
  {
    std::array<double_v_t,FunctionSpaceType::nDofsPerElement()> spatiallyVaryingPrefactorElementalValues;
    spatiallyVaryingPrefactor_->getElementValues(elementNoLocal, spatiallyVaryingPrefactorElementalValues);
    spatiallyVaryingPrefactor = functionSpace->interpolateValueInElement(spatiallyVaryingPrefactorElementalValues, xi);

    diffusionTensor = diffusionTensor * spatiallyVaryingPrefactor;
    if (VLOG_IS_ON(1))
      VLOG(1) << "elementNoLocal " << elementNoLocal << ", factor: " << spatiallyVaryingPrefactor << ", scaled diffusionTensor: " << diffusionTensor;
  }

  // check if diffusion tensor contains nan values
  if (D == 3)
  {
    if (!MathUtility::isFinite(diffusionTensor[0]) || !MathUtility::isFinite(diffusionTensor[1]) || !MathUtility::isFinite(diffusionTensor[2])
      || !MathUtility::isFinite(diffusionTensor[3]) || !MathUtility::isFinite(diffusionTensor[4]) || !MathUtility::isFinite(diffusionTensor[5])
      || !MathUtility::isFinite(diffusionTensor[6]) || !MathUtility::isFinite(diffusionTensor[7]) || !MathUtility::isFinite(diffusionTensor[8]))
    {
      LOG(ERROR) << "Directional diffusion tensor contains nans or infs.";
      LOG(INFO) << "elementNoLocal: " << elementNoLocal << ", xi: " << xi;
      LOG(INFO) << "directionVector: " << directionVector << " elemental direction values: " << elementalValues;
      LOG(INFO) << ", spatiallyVaryingPrefactor: " << spatiallyVaryingPrefactor;
      LOG(INFO) << "diffusionTensor from settings: " << std::endl << this->diffusionTensor_.value(elementNoLocal);
      LOG(INFO) << "additionalDiffusionTensor: " << std::endl << this->additionalDiffusionTensor_.value(elementNoLocal);
      LOG(INFO) << "resulting diffusion tensor in direction " << directionVector << ":" << std::endl << diffusionTensor;
     }
  }
  else if (D == 2)
  {
    if (!MathUtility::isFinite(diffusionTensor[0]) || !MathUtility::isFinite(diffusionTensor[1]) || !MathUtility::isFinite(diffusionTensor[2])
      || !MathUtility::isFinite(diffusionTensor[3]))
    {
      LOG(ERROR) << "Directional diffusion tensor contains nans or infs.";
      LOG(INFO) << "elementNoLocal: " << elementNoLocal << ", xi: " << xi;
      LOG(INFO) << "directionVector: " << directionVector << " elemental direction values: " << elementalValues;
      LOG(INFO) << "diffusionTensor: " << diffusionTensor << ", this->additionalDiffusionTensor: " << this->additionalDiffusionTensor_.value(elementNoLocal)
        << ", spatiallyVaryingPrefactor: " << spatiallyVaryingPrefactor;
     }
  }
  else if (D == 1)
  {
    if (!MathUtility::isFinite(diffusionTensor[0]))
    {
      LOG(ERROR) << "Directional diffusion tensor contains nans.";
      LOG(INFO) << "elementNoLocal: " << elementNoLocal << ", xi: " << xi;
      LOG(INFO) << "directionVector: " << directionVector << " elemental direction values: " << elementalValues;
      LOG(INFO) << "diffusionTensor: " << diffusionTensor << ", this->additionalDiffusionTensor: " << this->additionalDiffusionTensor_.value(elementNoLocal)
        << ", spatiallyVaryingPrefactor: " << spatiallyVaryingPrefactor;
     }
  }

  return diffusionTensor;
}

}  // namespace
