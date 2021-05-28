#include "mesh/mapping_between_meshes/mapping/01_implementation.h"

#include "control/diagnostic_tool/performance_measurement.h"

#include "utility/vector_operators.h"
#include "control/dihu_context.h"
#include "mesh/type_traits.h"
#include "mesh/mapping_between_meshes/manager/04_manager.h"
#include "mesh/mapping_between_meshes/manager/target_element_no_estimator.h"

//#define OUTPUT_INTERPOLATION_LEAP    // debugging output

namespace MappingBetweenMeshes
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  template<int nComponentsSource, int nComponentsTarget>
void MappingBetweenMeshesImplementation<FunctionSpaceSourceType, FunctionSpaceTargetType>::
mapLowToHighDimension(
  FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableSource, int componentNoSource,
  FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponentsTarget> &fieldVariableTarget, int componentNoTarget,
  FieldVariable::FieldVariable<FunctionSpaceTargetType,1> &targetFactorSum)
{
  assert(componentNoSource >= 0 && componentNoSource < nComponentsSource);
  assert(componentNoTarget >= 0 && componentNoTarget < nComponentsTarget);

  const dof_no_t nDofsLocalSource = fieldVariableSource.functionSpace()->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  std::vector<double> sourceValues;
  fieldVariableSource.getValuesWithoutGhosts(componentNoSource, sourceValues);

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "map " << fieldVariableSource.name() << "." << componentNoSource <<
      " (" << fieldVariableSource.functionSpace()->meshName() << ") -> " << fieldVariableTarget.name() << "." << componentNoTarget
      << " (" << fieldVariableTarget.functionSpace()->meshName() << ")";

    VLOG(1) << "source has " << nDofsLocalSource << " local dofs";
    VLOG(1) << fieldVariableSource;
    VLOG(1) << "extracted source values: " << sourceValues;
  }

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
    // if source dof is outside of target mesh, do nothing
    if (!this->targetMappingInfo_[sourceDofNoLocal].mapThisDof)
      continue;

    // get value
    double sourceValue = sourceValues[sourceDofNoLocal];

    // loop over target elements that will be affected by this source value
    for (int targetElementIndex = 0; targetElementIndex < this->targetMappingInfo_[sourceDofNoLocal].targetElements.size(); targetElementIndex++)
    {
      const typename MappingBetweenMeshesConstruct<FunctionSpaceSourceType,FunctionSpaceTargetType>::targetDof_t::element_t &targetElement
        = this->targetMappingInfo_[sourceDofNoLocal].targetElements[targetElementIndex];

      // store the value to the target function space
      element_no_t targetElementNoLocal = targetElement.elementNoLocal;
      std::array<double,nDofsPerTargetElement> targetValues = targetElement.scalingFactors * sourceValue;

      // determine dof nos of target element
      std::array<dof_no_t,nDofsPerTargetElement> dofNosLocal;
      for (int dofIndex = 0; dofIndex < nDofsPerTargetElement; dofIndex++)
      {
        dofNosLocal[dofIndex] = fieldVariableTarget.functionSpace()->getDofNo(targetElementNoLocal, dofIndex);
      }

      fieldVariableTarget.template setValues<nDofsPerTargetElement>(componentNoTarget, dofNosLocal, targetValues, ADD_VALUES);
      targetFactorSum.template setValues<nDofsPerTargetElement>(dofNosLocal, targetElement.scalingFactors, ADD_VALUES);

      if (VLOG_IS_ON(2))
      {
        VLOG(2) << "  source dof " << sourceDofNoLocal << ", value: " << sourceValue << ", scaling factors: " << targetElement.scalingFactors
          << ", targetValues: " << targetValues
          << ", targetElementNoLocal: " << targetElementNoLocal << ", target dofs: " << dofNosLocal;
      }
    }
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  template<int nComponents>
void MappingBetweenMeshesImplementation<FunctionSpaceSourceType, FunctionSpaceTargetType>::
mapLowToHighDimension(
  FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponents> &fieldVariableSource,
  FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponents> &fieldVariableTarget,
  FieldVariable::FieldVariable<FunctionSpaceTargetType,1> &targetFactorSum
)
{
  const dof_no_t nDofsLocalSource = fieldVariableSource.functionSpace()->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  std::vector<VecD<nComponents>> sourceValues;
  fieldVariableSource.getValuesWithoutGhosts(sourceValues);

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "map " << fieldVariableSource.name() << " (" << fieldVariableSource.functionSpace()->meshName()
      << ") -> " << fieldVariableTarget.name() << " (" << fieldVariableTarget.functionSpace()->meshName() << ")";

    VLOG(1) << "source has " << nDofsLocalSource << " local dofs";
    VLOG(1) << fieldVariableSource;
    VLOG(1) << "extracted source values: " << sourceValues;
  }

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
    VLOG(1) << "source dof " << sourceDofNoLocal << " map: " << this->targetMappingInfo_[sourceDofNoLocal].mapThisDof
      << ", has " << this->targetMappingInfo_[sourceDofNoLocal].targetElements.size() << " entry in target elements";

    // if target dof is outside of source mesh
    if (!this->targetMappingInfo_[sourceDofNoLocal].mapThisDof)
      continue;

    // get value
    VecD<nComponents> sourceValue = sourceValues[sourceDofNoLocal];

    // loop over target elements that will be affected by this source value
    for (int targetElementIndex = 0; targetElementIndex < this->targetMappingInfo_[sourceDofNoLocal].targetElements.size(); targetElementIndex++)
    {
      const typename MappingBetweenMeshesConstruct<FunctionSpaceSourceType,FunctionSpaceTargetType>::targetDof_t::element_t &targetElement
        = this->targetMappingInfo_[sourceDofNoLocal].targetElements[targetElementIndex];

      // store the value to the target function space
      element_no_t targetElementNoLocal = targetElement.elementNoLocal;
      std::array<VecD<nComponents>,nDofsPerTargetElement> targetValues = targetElement.scalingFactors * sourceValue;

      // determine dof nos of target element
      std::array<dof_no_t,nDofsPerTargetElement> dofNosLocal;
      for (int dofIndex = 0; dofIndex < nDofsPerTargetElement; dofIndex++)
      {
        dofNosLocal[dofIndex] = fieldVariableTarget.functionSpace()->getDofNo(targetElementNoLocal, dofIndex);

        if (dofNosLocal[dofIndex] >= fieldVariableTarget.functionSpace()->nDofsLocalWithGhosts() 
            || dofNosLocal[dofIndex] >= targetFactorSum.functionSpace()->nDofsLocalWithGhosts())
        {
          LOG(FATAL) << "Dof no " << dofNosLocal[dofIndex] << " out of range, \"" << fieldVariableTarget.functionSpace()->meshName() << "\" has "
            << fieldVariableTarget.functionSpace()->nDofsLocalWithGhosts() << " local dofs with ghosts, \"" << targetFactorSum.functionSpace()->meshName() << "\" has "
            << targetFactorSum.functionSpace()->nDofsLocalWithGhosts() << " local dofs with ghosts. Mapping "
            << fieldVariableSource.name() << " (" << fieldVariableSource.functionSpace()->meshName()
            << ") -> " << fieldVariableTarget.name() << " (" << fieldVariableTarget.functionSpace()->meshName() << "), targetElementNoLocal: " 
            << targetElementNoLocal << "/" << fieldVariableTarget.functionSpace()->nElementsLocal() << ", dofIndex: " << dofIndex << "/" << nDofsPerTargetElement;
        }
      }

      fieldVariableTarget.template setValues<nDofsPerTargetElement>(dofNosLocal, targetValues, ADD_VALUES);
      targetFactorSum.template setValues<nDofsPerTargetElement>(dofNosLocal, targetElement.scalingFactors, ADD_VALUES);

      if (VLOG_IS_ON(2))
      {
        VLOG(2) << "  source dof " << sourceDofNoLocal << ", value: " << sourceValue << ", scaling factors: " << targetElement.scalingFactors
          << ", targetValues: " << targetValues
          << ", targetElementNoLocal: " << targetElementNoLocal << ", target dofs: " << dofNosLocal;
      }
    }
  }
}

//! map data between all components of the field variables in the source and target function spaces
//! note the swapped template arguments of MappingBetweenMeshesImplementation, in this whole method the mapping is still source->target,
//! but the function space names are different than in the other methods.
template<typename FunctionSpaceTargetType, typename FunctionSpaceSourceType>
template<int nComponents>
void MappingBetweenMeshesImplementation<FunctionSpaceTargetType, FunctionSpaceSourceType>::
mapHighToLowDimension(
  FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponents> &fieldVariableSource,
  FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponents> &fieldVariableTarget
)
{
  const dof_no_t nDofsLocalTarget = fieldVariableTarget.functionSpace()->nDofsLocalWithoutGhosts();
  const int nDofsPerSourceElement = FunctionSpaceSourceType::nDofsPerElement();

  if (VLOG_IS_ON(1))
  {
    std::vector<VecD<nComponents>> sourceValues;
    fieldVariableSource.getValuesWithoutGhosts(sourceValues);

    VLOG(1) << "map " << fieldVariableSource.name() << " (" << fieldVariableSource.functionSpace()->meshName()
      << ") -> " << fieldVariableTarget.name() << " (" << fieldVariableTarget.functionSpace()->meshName() << ")";

    VLOG(1) << "source has " << nDofsLocalTarget << " local dofs";
    VLOG(1) << fieldVariableSource;
    VLOG(1) << "extracted source values: " << sourceValues;
  }

  fieldVariableSource.zeroGhostBuffer();
  fieldVariableSource.finishGhostManipulation();
  fieldVariableSource.startGhostManipulation();

#ifdef OUTPUT_INTERPOLATION_LEAP
  VecD<nComponents> previousTargetValue;
  element_no_t previousSourceElementNoLocal;
  std::array<double,FunctionSpaceSourceType::nDofsPerElement()> previousScalingFactors;
#endif

  LOG(DEBUG) << "mapHighToLowDimension " << fieldVariableSource.name() << " (" << fieldVariableSource.functionSpace()->meshName()
      << ") -> " << fieldVariableTarget.name() << " (" << fieldVariableTarget.functionSpace()->meshName() << ")";

  // this mapping direction corresponds to simple interpolation in the source mesh

  // visualization for 1D-1D: s=source, t=target
  // s--t--------s-----t-----s

  // loop over all local dofs of the target functionSpace, which was the source function space when the mapping was initialized
  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    // if source dof is outside of target mesh
    if (!this->targetMappingInfo_[targetDofNoLocal].mapThisDof)
      continue;

    // the first set of surrounding nodes (targetElements[0]) is enough
    const typename MappingBetweenMeshesConstruct<FunctionSpaceTargetType,FunctionSpaceSourceType>::targetDof_t::element_t &targetElement
      = this->targetMappingInfo_[targetDofNoLocal].targetElements[0];

    element_no_t sourceElementNoLocal = targetElement.elementNoLocal;

    // get source values of the element where targetDofNoLocal is in
    std::array<VecD<nComponents>,nDofsPerSourceElement> sourceValues;
    fieldVariableSource.getElementValues(sourceElementNoLocal, sourceValues);

    double scalingFactorsSum = 0;
    for (int i = 0; i < nDofsPerSourceElement; i++)
    {
      scalingFactorsSum += targetElement.scalingFactors[i];
    }

    if (fabs(scalingFactorsSum-1.0) > 1e-10)
      LOG(ERROR) << "Scaling factors do not sum to 1, scalingFactorsSum: " << scalingFactorsSum << ", scalingFactors: " << targetElement.scalingFactors;

    VecD<nComponents> targetValue = sourceValues * targetElement.scalingFactors;
    fieldVariableTarget.setValue(targetDofNoLocal, targetValue, INSERT_VALUES);

    LOG(DEBUG) << fieldVariableTarget.functionSpace()->meshName() << " dof " << targetDofNoLocal << " value = " << targetValue << " = " << sourceValues << " * " << targetElement.scalingFactors << ", interpolation in element " << sourceElementNoLocal << " of " << fieldVariableSource.functionSpace()->meshName();

#ifdef OUTPUT_INTERPOLATION_LEAP
    if (targetDofNoLocal > 0)
    {
      double differenceToPrevious = MathUtility::distance<nComponents>(targetValue, previousTargetValue);
      LOG(INFO) << "  target dof " << targetDofNoLocal << ", source element no local: " << sourceElementNoLocal << " differenceToPrevious=" << differenceToPrevious;
      if (differenceToPrevious > 1)
      {
        LOG(WARNING) << "In mapping " << fieldVariableSource.name() << " (" << fieldVariableSource.functionSpace()->meshName()
          << ") -> " << fieldVariableTarget.name() << " (" << fieldVariableTarget.functionSpace()->meshName() << "), differenceToPrevious=" << differenceToPrevious;
        LOG(WARNING) << "  target dof " << targetDofNoLocal << ", source element no local: " << sourceElementNoLocal
        << ", " << sourceValues.size() << " sourceValues: " << sourceValues
        << ",\n scaling factors: " << targetElement.scalingFactors << " (scalingFactorsSum=" << scalingFactorsSum << "),\n previous target value: " << previousTargetValue
        << ", target value: " << targetValue << ", previous source element: " << previousSourceElementNoLocal
        << "\npreviousScalingFactors: " << previousScalingFactors;
      }
    }
    if (VLOG_IS_ON(2))
    {
      VLOG(2)
        << "  target dof " << targetDofNoLocal << ", source element no local: " << sourceElementNoLocal
        << ", sourceValues: " << sourceValues
        << ", scaling factors: " << targetElement.scalingFactors << " (scalingFactorsSum=" << scalingFactorsSum << "), target value: " << targetValue;
    }
    // store quantities from last dof for debugging output
    previousTargetValue = targetValue;
    previousSourceElementNoLocal = sourceElementNoLocal;
    previousScalingFactors = targetElement.scalingFactors;
#endif
  }
}

//! map data between specific components of the field variables in the source and target function spaces
//! note the swapped template arguments of MappingBetweenMeshesImplementation, in this whole method the mapping is still source->target,
//! but the function space names are different than in the other methods.
template<typename FunctionSpaceTargetType, typename FunctionSpaceSourceType>
template<int nComponentsSource, int nComponentsTarget>
void MappingBetweenMeshesImplementation<FunctionSpaceTargetType, FunctionSpaceSourceType>::
mapHighToLowDimension(
  FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableSource, int componentNoSource,
  FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponentsTarget> &fieldVariableTarget, int componentNoTarget
)
{
  assert(componentNoSource >= 0 && componentNoSource < nComponentsSource);
  assert(componentNoTarget >= 0 && componentNoTarget < nComponentsTarget);

  const dof_no_t nDofsLocalTarget = fieldVariableTarget.functionSpace()->nDofsLocalWithoutGhosts();
  const int nDofsPerSourceElement = FunctionSpaceSourceType::nDofsPerElement();

  if (VLOG_IS_ON(1))
  {
    std::vector<double> sourceValues;
    fieldVariableSource.getValuesWithoutGhosts(componentNoSource, sourceValues);

    VLOG(1) << "map " << fieldVariableSource.name() << "." << componentNoSource
      << " (" << fieldVariableSource.functionSpace()->meshName()
      << ") -> " << fieldVariableTarget.name() << "." << componentNoTarget
      << " (" << fieldVariableTarget.functionSpace()->meshName() << ")";

    VLOG(1) << "source has " << nDofsLocalTarget << " local dofs";
    VLOG(1) << fieldVariableSource;
    VLOG(1) << "extracted source values: " << sourceValues;
  }

  // visualization for 1D-1D: s=source, t=target
  // s--t--------s-----t-----s

  // loop over all local dofs of the target functionSpace, which was the source function space when the mapping was initialized
  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    // if source dof is outside of target mesh
    if (!this->targetMappingInfo_[targetDofNoLocal].mapThisDof)
      continue;

    // the first set of surrounding nodes (targetElements[0]) is enough
    const typename MappingBetweenMeshesConstruct<FunctionSpaceTargetType,FunctionSpaceSourceType>::targetDof_t::element_t &targetElement
      = this->targetMappingInfo_[targetDofNoLocal].targetElements[0];

    element_no_t sourceElementNoLocal = targetElement.elementNoLocal;

    double scalingFactorsSum = 0;
    for (int i = 0; i < nDofsPerSourceElement; i++)
    {
      scalingFactorsSum += targetElement.scalingFactors[i];
    }

    // get source values of the element where targetDofNoLocal is in
    std::array<double,nDofsPerSourceElement> sourceValues;
    fieldVariableSource.getElementValues(componentNoSource, sourceElementNoLocal, sourceValues);

    double targetValue = 0;
    for (int i = 0; i < nDofsPerSourceElement; i++)
    {
      targetValue += sourceValues[i] * targetElement.scalingFactors[i] / scalingFactorsSum;
    }
    fieldVariableTarget.setValue(componentNoTarget, targetDofNoLocal, targetValue, INSERT_VALUES);

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "  target dof " << targetDofNoLocal << ", source element no local: " << sourceElementNoLocal
        << ", sourceValues: " << sourceValues
        << ", scaling factors: " << targetElement.scalingFactors << ", target value: " << targetValue;
    }
  }
}


}  // namespace
