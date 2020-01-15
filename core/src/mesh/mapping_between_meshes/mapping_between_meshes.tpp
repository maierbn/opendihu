#include "mesh/mapping_between_meshes/mapping_between_meshes.h"

#include "control/performance_measurement.h"

#include "utility/vector_operators.h"

namespace Mesh
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::
MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                     std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget,
                     double xiTolerance) :
  functionSpaceSource_(functionSpaceSource),
  functionSpaceTarget_(functionSpaceTarget),
  maxAllowedXiTolerance_(xiTolerance)
{
  // create the mapping

  Control::PerformanceMeasurement::start("durationComputeMappingBetweenMeshes");

  const dof_no_t nDofsLocalSource = functionSpaceSource->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  element_no_t elementNo = 0;
  int ghostMeshNo = 0;
  std::array<double,FunctionSpaceTargetType::dim()> xi;

  targetMappingInfo_.resize(nDofsLocalSource);

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();
    VLOG(1) << "geometry: " << functionSpaceSource->geometryField();
  }

  //VLOG(1) << "target meshPartition: " << *functionSpaceTarget->meshPartition();
  //VLOG(1) << "geometryField: " << functionSpaceTarget->geometryField();

  double xiToleranceBase = 1e-2;    // internal tolerance is 1e-3

  bool mappingSucceeded = true;
  bool startSearchInCurrentElement = false;
  int nSourceDofsOutsideTargetMesh = 0;

  // visualization for 1D-1D: s=source, t=target
  // t--s--------t-----s-----t

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
    // determine information how to map a source value to the target mesh
    targetDof_t targetMappingInfo;

    // get node position of the source dof

    //dof_no_t sourceDofNoGlobal = functionSpaceTarget->meshPartition()->getDofNoGlobalPetsc(sourceDofNoLocal);
    Vec3 position = functionSpaceSource->getGeometry(sourceDofNoLocal);

    double xiTolerance = xiToleranceBase;
    int nTries = 0;
    int nTriesMax = 3;
    int startElementNo = elementNo;
    for(nTries = 0; nTries < nTriesMax; nTries++)
    {
      // find element no in the target mesh where the position is
      elementNo = startElementNo;
      if (functionSpaceTarget->findPosition(position, elementNo, ghostMeshNo, xi, startSearchInCurrentElement, xiTolerance))
      {
        // If there was no prescribed maxAllowedXiTolerance_, set the new xiTolerance to the value from which the current search succeeded,
        // because it is assumed the the current two meshes are located to each other so mismatching that this tolerance is enough.
        if (maxAllowedXiTolerance_ == 0)
        {
          xiToleranceBase = xiTolerance;
        }

        targetMappingInfo.mapThisDof = true;
        break;
      }
      else
      {
        xiTolerance *= 2;
        if (maxAllowedXiTolerance_ != 0 && xiTolerance > maxAllowedXiTolerance_)
        {
          break;
        }
        LOG(DEBUG) << "Try again with xiTolerance = " << xiTolerance;
      }
    }

    if (nTries == nTriesMax || (maxAllowedXiTolerance_ != 0 && xiTolerance > maxAllowedXiTolerance_))
    {
      LOG(DEBUG) << "In mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
        << functionSpaceTarget->meshName() << "\", source dof local " << sourceDofNoLocal
        << " of mesh \"" << functionSpaceSource->meshName() << "\" at position " << position << " is outside of target mesh \""
        << functionSpaceTarget->meshName() << "\" with tolerance " << xiTolerance << ".";

      nSourceDofsOutsideTargetMesh++;
      targetMappingInfo.mapThisDof = false;
    }

    // store element no
    targetMappingInfo.elementNoLocal = elementNo;

    // determine factors how to distribute the value to the dofs of the target element

    // note: geometry value = sum over dofs of geometryValue_dof * phi_dof(xi)
    for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
    {
      VLOG(3) << "   phi_" << targetDofIndex << "(" << xi << ")=" << functionSpaceTarget->phi(targetDofIndex, xi);
      double phiContribution = functionSpaceTarget->phi(targetDofIndex, xi);

      if (fabs(phiContribution) < 1e-7)
      {
        if (phiContribution >= 0)
        {
          phiContribution = 1e-7;
        }
        else
        {
          phiContribution = -1e-7;
        }
      }
      targetMappingInfo.scalingFactors[targetDofIndex] = phiContribution;
    }

    targetMappingInfo_[sourceDofNoLocal] = targetMappingInfo;

    if (VLOG_IS_ON(2))
    {
      double scalingFactorsSum = 0;
      for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
      {
        scalingFactorsSum += targetMappingInfo_[sourceDofNoLocal].scalingFactors[targetDofIndex];
      }
      VLOG(3) << "  source dof local " << sourceDofNoLocal << ", pos: " << position << ", xi: " << xi
        << ", element no: " << targetMappingInfo.elementNoLocal << ", scaling factors: " << targetMappingInfo_[sourceDofNoLocal].scalingFactors
        << ", sum: " << scalingFactorsSum;
    }

    // next time when searching for the target element, start search from previous element
    startSearchInCurrentElement = true;
  }

  Control::PerformanceMeasurement::stop("durationComputeMappingBetweenMeshes");

  if (!mappingSucceeded)
  {
    LOG(ERROR) << "Could not create mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
      << functionSpaceTarget->meshName() << "\".";
    LOG(FATAL) << "end";
  }
  else
  {
    LOG(DEBUG) << "Successfully initialized mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
      << functionSpaceTarget->meshName() << "\", " << nSourceDofsOutsideTargetMesh << "/" << nDofsLocalSource << " source dofs are outside the target mesh.";
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  template<int nComponentsSource, int nComponentsTarget>
void MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::
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
    // if source dof is outside of target mesh
    if (!targetMappingInfo_[sourceDofNoLocal].mapThisDof)
      continue;

    // get value
    double sourceValue = sourceValues[sourceDofNoLocal];

    // store the value to the target function space
    element_no_t targetElementNoLocal = targetMappingInfo_[sourceDofNoLocal].elementNoLocal;
    std::array<double,nDofsPerTargetElement> targetValues = targetMappingInfo_[sourceDofNoLocal].scalingFactors * sourceValue;

    // determine dof nos of target element
    std::array<dof_no_t,nDofsPerTargetElement> dofNosLocal;
    for (int dofIndex = 0; dofIndex < nDofsPerTargetElement; dofIndex++)
    {
      dofNosLocal[dofIndex] = fieldVariableTarget.functionSpace()->getDofNo(targetElementNoLocal, dofIndex);
    }

    fieldVariableTarget.template setValues<nDofsPerTargetElement>(componentNoTarget, dofNosLocal, targetValues, ADD_VALUES);
    targetFactorSum.template setValues<nDofsPerTargetElement>(dofNosLocal, targetMappingInfo_[sourceDofNoLocal].scalingFactors, ADD_VALUES);

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "  source dof " << sourceDofNoLocal << ", value: " << sourceValue << ", scaling factors: " << targetMappingInfo_[sourceDofNoLocal].scalingFactors
        << ", targetValues: " << targetValues
        << ", targetElementNoLocal: " << targetElementNoLocal << ", target dofs: " << dofNosLocal;
    }
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  template<int nComponents>
void MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::
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
    // if target dof is outside of source mesh
    if (!targetMappingInfo_[sourceDofNoLocal].mapThisDof)
      continue;

    // get value
    VecD<nComponents> sourceValue = sourceValues[sourceDofNoLocal];

    // store the value to the target function space
    element_no_t targetElementNoLocal = targetMappingInfo_[sourceDofNoLocal].elementNoLocal;
    std::array<VecD<nComponents>,nDofsPerTargetElement> targetValues = targetMappingInfo_[sourceDofNoLocal].scalingFactors * sourceValue;

    // determine dof nos of target element
    std::array<dof_no_t,nDofsPerTargetElement> dofNosLocal;
    for (int dofIndex = 0; dofIndex < nDofsPerTargetElement; dofIndex++)
    {
      dofNosLocal[dofIndex] = fieldVariableTarget.functionSpace()->getDofNo(targetElementNoLocal, dofIndex);
    }

    fieldVariableTarget.template setValues<nDofsPerTargetElement>(dofNosLocal, targetValues, ADD_VALUES);
    targetFactorSum.template setValues<nDofsPerTargetElement>(dofNosLocal, targetMappingInfo_[sourceDofNoLocal].scalingFactors, ADD_VALUES);

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "  source dof " << sourceDofNoLocal << ", value: " << sourceValue << ", scaling factors: " << targetMappingInfo_[sourceDofNoLocal].scalingFactors
        << ", targetValues: " << targetValues
        << ", targetElementNoLocal: " << targetElementNoLocal << ", target dofs: " << dofNosLocal;
    }
  }
}

//! map data between all components of the field variables in the source and target function spaces
//! note the swapped template arguments of MappingBetweenMeshes, in this whole method the mapping is still source->target,
//! but the function space names are different than in the other methods.
template<typename FunctionSpaceTargetType, typename FunctionSpaceSourceType>
template<int nComponents>
void MappingBetweenMeshes<FunctionSpaceTargetType, FunctionSpaceSourceType>::
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

  // visualization for 1D-1D: s=source, t=target
  // s--t--------s-----t-----s

  // loop over all local dofs of the target functionSpace, which was the source function space when the mapping was initialized
  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    // if source dof is outside of target mesh
    if (!targetMappingInfo_[targetDofNoLocal].mapThisDof)
      continue;

    element_no_t sourceElementNoLocal = targetMappingInfo_[targetDofNoLocal].elementNoLocal;

    // get source values of the element where targetDofNoLocal is in
    std::array<VecD<nComponents>,nDofsPerSourceElement> sourceValues;
    fieldVariableSource.getElementValues(sourceElementNoLocal, sourceValues);

    VecD<nComponents> targetValue = sourceValues * targetMappingInfo_[targetDofNoLocal].scalingFactors;
    fieldVariableTarget.setValue(targetDofNoLocal, targetValue);

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "  target dof " << targetDofNoLocal << ", source element no local: " << sourceElementNoLocal
        << ", sourceValues: " << sourceValues
        << ", scaling factors: " << targetMappingInfo_[targetDofNoLocal].scalingFactors << ", target value: " << targetValue;
    }
  }
}

//! map data between specific components of the field variables in the source and target function spaces
//! note the swapped template arguments of MappingBetweenMeshes, in this whole method the mapping is still source->target,
//! but the function space names are different than in the other methods.
template<typename FunctionSpaceTargetType, typename FunctionSpaceSourceType>
template<int nComponentsSource, int nComponentsTarget>
void MappingBetweenMeshes<FunctionSpaceTargetType, FunctionSpaceSourceType>::
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
    if (!targetMappingInfo_[targetDofNoLocal].mapThisDof)
      continue;

    element_no_t sourceElementNoLocal = targetMappingInfo_[targetDofNoLocal].elementNoLocal;

    // get source values of the element where targetDofNoLocal is in
    std::array<double,nDofsPerSourceElement> sourceValues;
    fieldVariableSource.getElementValues(componentNoSource, sourceElementNoLocal, sourceValues);

    double targetValue = 0;
    for (int i = 0; i < nDofsPerSourceElement; i++)
    {
      targetValue += sourceValues[i] * targetMappingInfo_[targetDofNoLocal].scalingFactors[i];
    }
    fieldVariableTarget.setValue(componentNoTarget, targetDofNoLocal, targetValue, INSERT_VALUES);

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "  target dof " << targetDofNoLocal << ", source element no local: " << sourceElementNoLocal
        << ", sourceValues: " << sourceValues
        << ", scaling factors: " << targetMappingInfo_[targetDofNoLocal].scalingFactors << ", target value: " << targetValue;
    }
  }
}

}  // namespace
