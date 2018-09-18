#include "data_management/diffusion_tensor_field_variable.h"

namespace Data
{


template<typename FunctionSpaceType>
void DiffusionTensorFieldVariable<FunctionSpaceType>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction)
{
  LOG(DEBUG) << "initialize Diffusion tensor using field variable";
  direction_ = direction;
}

template<typename FunctionSpaceType>
MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> DiffusionTensor<FunctionSpaceType>::
diffusionTensor(element_no_t elementNoLocal, const std::array<double,FunctionSpaceType::dim()> xi) const
{
  const int D = FunctionSpaceType::dim();

  // get the values at dofs of the element
  std::array<Vec3,FunctionSpaceType::nDofsPerElement()> elementalValues;
  direction_->getElementValues(elementNo, elementalValues);

  // get the function space from field variable
  std::shared_ptr<FunctionSpaceType> functionSpace = direction_->functionSpace();

  // get the interpolated value at xi coordinates inside the element
  Vec3 directionVector = functionSpace->template interpolateValueInElement<3>(elementalValues,xi) const;

  MathUtility::Matrix<D,D> diffusionTensor;

  if (FunctionSpaceType::dim() == 1)
  {
    diffusionTensor = {1.0};
  }
  else if (FunctionSpaceType::dim() == 2)
  {
    diffusionTensorExtraCellular = {1.0, 0.0,
                                    0.0, 0.0};
  }
  else if (FunctionSpaceType::dim() == 3)
  {
    diffusionTensorExtraCellular = {1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0};
  }

  MathUtility::rotateMatrix(diffusionTensorExtraCellular, directionVector);

  if (multidomainNCompartments_ > 0)
  {
    MathUtility::Matrix<D,D> diffusionTensorExtraCellular;  // sigma_e

    if (FunctionSpaceType::dim() == 1)
    {
      diffusionTensorExtraCellular = {1.0};
    }
    else if (FunctionSpaceType::dim() == 2)
    {
      diffusionTensorExtraCellular = {2.0, 0.0,
                                      0.0, 1.0};
    }
    else if (FunctionSpaceType::dim() == 3)
    {
      diffusionTensorExtraCellular = {1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 2.0};
    }

    MathUtility::rotateMatrix(diffusionTensorExtraCellular, directionVector);

    // add extracellular diffusion tensor
    diffusionTensor = diffusionTensor + diffusionTensorExtraCellular;
  }

  return diffusionTensor;
}

}  // namespace
