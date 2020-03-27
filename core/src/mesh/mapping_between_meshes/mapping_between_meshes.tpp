#include "mesh/mapping_between_meshes/mapping_between_meshes.h"

#include "control/diagnostic_tool/performance_measurement.h"

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
  const dof_no_t nDofsLocalTarget = functionSpaceTarget->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  LOG(DEBUG) << "create MappingBetweenMeshes \"" << functionSpaceSource->meshName() << "\" and \""
        << functionSpaceTarget->meshName() << "\" xiTolerance: " << xiTolerance;

  element_no_t elementNo = 0;
  int ghostMeshNo = 0;
  std::array<double,FunctionSpaceTargetType::dim()> xi;

  targetMappingInfo_.resize(nDofsLocalSource);

  std::vector<bool> targetDofIsMappedTo(nDofsLocalTarget, false);   //< for every target dof if it will get a value from any source dof

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();
    VLOG(1) << "geometry: " << functionSpaceSource->geometryField();
  }

  //VLOG(1) << "target meshPartition: " << *functionSpaceTarget->meshPartition();
  //VLOG(1) << "geometryField: " << functionSpaceTarget->geometryField();

  if (xiTolerance <= 0)
    xiTolerance = 1e-2;
    
  bool startSearchInCurrentElement = false;
  int nSourceDofsOutsideTargetMesh = 0;
  double residual;

  // visualization for 1D-1D: s=source, t=target
  // t--s--------t-----s-----t

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
    // determine information how to map a source value to the target mesh
    targetDof_t targetMappingInfo;
    targetMappingInfo.targetElements.resize(1);

    // get node position of the source dof

    //dof_no_t sourceDofNoGlobal = functionSpaceTarget->meshPartition()->getDofNoGlobalPetsc(sourceDofNoLocal);
    Vec3 position = functionSpaceSource->getGeometry(sourceDofNoLocal);

    // find element no in the target mesh where the position is
    if (functionSpaceTarget->findPosition(position, elementNo, ghostMeshNo, xi, startSearchInCurrentElement, residual, xiTolerance))
    {
      targetMappingInfo.mapThisDof = true;
      VLOG(1) << "found at xi=" << xi << ", elementNo: " << elementNo << ", xiTolerance=" << xiTolerance << ", residual: " << residual;
    }
    else
    {
      LOG(DEBUG) << "In mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
        << functionSpaceTarget->meshName() << "\", source dof local " << sourceDofNoLocal
        << " of mesh \"" << functionSpaceSource->meshName() << "\" at position " << position << " is outside of target mesh \""
        << functionSpaceTarget->meshName() << "\" with tolerance " << xiTolerance << ". Try increasing parameter \"xiTolerance\".";

      nSourceDofsOutsideTargetMesh++;
      targetMappingInfo.mapThisDof = false;
    }

    // store element no
    targetMappingInfo.targetElements[0].elementNoLocal = elementNo;

    std::array<dof_no_t,FunctionSpaceTargetType::nDofsPerElement()> targetDofNos = functionSpaceTarget->getElementDofNosLocal(elementNo);

    // determine factors how to distribute the value to the dofs of the target element

    // note: geometry value = sum over dofs of geometryValue_dof * phi_dof(xi)
    for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
    {
      VLOG(3) << "   phi_" << targetDofIndex << "(" << xi << ")=" << functionSpaceTarget->phi(targetDofIndex, xi);
      double phiContribution = functionSpaceTarget->phi(targetDofIndex, xi);

      // if phi is close to zero, set to 1e-14, this is practically zero, but it is still possible to divide by it in case the dof does not get any other contribution
      if (fabs(phiContribution) < 1e-14)
      {
        if (phiContribution >= 0)
        {
          phiContribution = 1e-14;
        }
        else
        {
          phiContribution = -1e-14;
        }
      }
      else
      {
        dof_no_t targetDofNoLocal = targetDofNos[targetDofIndex];

        // if this dof is local, store information that this target dof will get a value in the mapping, i.e. there is a source dof that influences the mapped value of the target dof
        if (targetDofNoLocal < nDofsLocalTarget)
        {
          targetDofIsMappedTo[targetDofNoLocal] = true;
        }
      }

      targetMappingInfo.targetElements[0].scalingFactors[targetDofIndex] = phiContribution;
    }
    

    targetMappingInfo_[sourceDofNoLocal] = targetMappingInfo;

    if (VLOG_IS_ON(2))
    {
      double scalingFactorsSum = 0;
      for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
      {
        scalingFactorsSum += targetMappingInfo_[sourceDofNoLocal].targetElements[0].scalingFactors[targetDofIndex];
      }
      VLOG(2) << "  source dof local " << sourceDofNoLocal << ", pos: " << position << ", xi: " << xi
        << ", element no: " << targetMappingInfo.targetElements[0].elementNoLocal << ", scaling factors: " << targetMappingInfo_[sourceDofNoLocal].targetElements[0].scalingFactors
        << ", sum: " << scalingFactorsSum;
    }

    // next time when searching for the target element, start search from previous element (the current element), because here we found the point and the next point is probably in one of the neighbouring elements
    startSearchInCurrentElement = true;
  }

  // find target dofs that do not appear in any targetMappingInfo and therefore will so far not receive any value when mapping from source to target

  // count number of target dofs that have no source dof mapped
  int nTargetDofsNotMapped = 0;

  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal < nDofsLocalTarget; targetDofNoLocal++)
  {
    if (!targetDofIsMappedTo[targetDofNoLocal])
    {
      nTargetDofsNotMapped++;
    }
  }

  if (FunctionSpaceSourceType::dim() >= FunctionSpaceTargetType::dim() && nTargetDofsNotMapped > 0)
  {
    LOG(DEBUG) << nTargetDofsNotMapped << " target dofs have no source dofs that would contribute values. Source FunctionSpace dim: "
      << FunctionSpaceSourceType::dim() << " >= target FunctionSpace dim: " << FunctionSpaceTargetType::dim() << ". Now fixing.";

    std::set<dof_no_t> targetDofNoLocalNotFixed;    // collect all dofs that are still not fixed

    // loop over all elements in the target function space
    for (element_no_t targetElementNoLocal = 0; targetElementNoLocal < functionSpaceTarget->nElementsLocal(); targetElementNoLocal++)
    {
      std::array<dof_no_t,FunctionSpaceTargetType::nDofsPerElement()> targetDofNos = functionSpaceTarget->getElementDofNosLocal(targetElementNoLocal);

      // loop over dofs of target element
      for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
      {
        dof_no_t targetDofNoLocal = targetDofNos[targetDofIndex];

        // if this is a ghost dof, do not handle it
        if (targetDofNoLocal >= nDofsLocalTarget)
          continue;

        Vec3 position = functionSpaceTarget->getGeometry(targetDofNoLocal);
        LOG(DEBUG) << " e" << targetElementNoLocal << " i" << targetDofIndex << ", targetDofIsMappedTo[" << targetDofNoLocal << "]: " << targetDofIsMappedTo[targetDofNoLocal] << ", position: " << position;

        // if target dof is not being mapped to by any source dof, simply initiate interpolation of the source mesh to this dof
        if (!targetDofIsMappedTo[targetDofNoLocal])
        {
          // get one element of the target dof

          // find element and xi position in source mesh where target dof is located
          Vec3 position = functionSpaceTarget->getGeometry(targetDofNoLocal);
          element_no_t sourceElementNo = 0;
          bool startSearchInCurrentElement = false;
          std::array<double,FunctionSpaceSourceType::dim()> xiSource;

          //LOG(DEBUG) << "target (el." << targetElementNoLocal << ",index" << targetDofIndex << ") dof " << targetDofNoLocal << ", position: " << position
          //  << " is not mapped, now find element in source function space";

          if (functionSpaceSource->findPosition(position, sourceElementNo, ghostMeshNo, xiSource, startSearchInCurrentElement, residual, xiTolerance))
          {
            // get dofs of this source element
            std::array<dof_no_t,FunctionSpaceSourceType::nDofsPerElement()> sourceDofNos = functionSpaceSource->getElementDofNosLocal(sourceElementNo);

            LOG(DEBUG) << "at position " << position << " found source element " << sourceElementNo << ", xi " << xiSource << ", with dofs " << sourceDofNos;

            // loop over all the source dofs that will contribute to the value of the target dof
            for (int sourceDofIndex = 0; sourceDofIndex != FunctionSpaceSourceType::nDofsPerElement(); sourceDofIndex++)
            {
              dof_no_t sourceDofNoLocal = sourceDofNos[sourceDofIndex];

              // create new entry for the targetMappingInfo_[sourceDofNoLocal]
              typename targetDof_t::element_t targetElement;

              // set element no
              targetElement.elementNoLocal = targetElementNoLocal;

              // set scaling factors
              for (int i = 0; i < FunctionSpaceTargetType::nDofsPerElement(); i++)
              {
                targetElement.scalingFactors[i] = 0;
              }

              double phiContribution = functionSpaceSource->phi(sourceDofIndex, xiSource);
              targetElement.scalingFactors[targetDofIndex] = phiContribution;

              targetMappingInfo_[sourceDofNoLocal].targetElements.push_back(targetElement);

              //LOG(DEBUG) << "add scaling Factor " << phiContribution << " at targetDofIndex " << targetDofIndex << " of targetELement " << targetElementNoLocal
              //  << ", now, source dof " << sourceDofNoLocal << " has " << targetMappingInfo_[sourceDofNoLocal].targetElements.size() << " target elements.";

            }

            // now the target dof is fixed
            targetDofIsMappedTo[targetDofNoLocal] = true;
          }
          else
          {
            LOG(DEBUG) << "Could not find element of source dof for position " << position;
            targetDofNoLocalNotFixed.insert(targetDofNoLocal);
          }
        }
      }
    }
    LOG(DEBUG) << "after fixing target dofs by source mesh interpolation, " << targetDofNoLocalNotFixed.size()
      << " remaining targetDofNoLocalNotFixed: " << targetDofNoLocalNotFixed;
  }
  else
  {
    LOG(DEBUG) << nTargetDofsNotMapped << " target dofs have no source dofs that would contribute values. Source FunctionSpace dim: "
      << FunctionSpaceSourceType::dim() << ", target FunctionSpace dim: " << FunctionSpaceTargetType::dim() << ".";
  }

  Control::PerformanceMeasurement::stop("durationComputeMappingBetweenMeshes");

  if (nSourceDofsOutsideTargetMesh > 0)
  {
    LOG(INFO) << "Successfully initialized mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
      << functionSpaceTarget->meshName() << "\", " << nSourceDofsOutsideTargetMesh << "/" << nDofsLocalSource << " source dofs are outside the target mesh. "
      << "\"xiTolerance\": " << xiTolerance << ", total duration of all mappings: " << Control::PerformanceMeasurement::getDuration("durationComputeMappingBetweenMeshes") << " s";
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
    // if source dof is outside of target mesh, do nothing
    if (!targetMappingInfo_[sourceDofNoLocal].mapThisDof)
      continue;

    // get value
    double sourceValue = sourceValues[sourceDofNoLocal];

    // loop over target elements that will be affected by this source value
    for (int targetElementIndex = 0; targetElementIndex < targetMappingInfo_[sourceDofNoLocal].targetElements.size(); targetElementIndex++)
    {
      const typename targetDof_t::element_t &targetElement = targetMappingInfo_[sourceDofNoLocal].targetElements[targetElementIndex];

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

    // loop over target elements that will be affected by this source value
    for (int targetElementIndex = 0; targetElementIndex < targetMappingInfo_[sourceDofNoLocal].targetElements.size(); targetElementIndex++)
    {
      const typename targetDof_t::element_t &targetElement = targetMappingInfo_[sourceDofNoLocal].targetElements[targetElementIndex];

      // store the value to the target function space
      element_no_t targetElementNoLocal = targetElement.elementNoLocal;
      std::array<VecD<nComponents>,nDofsPerTargetElement> targetValues = targetElement.scalingFactors * sourceValue;

      // determine dof nos of target element
      std::array<dof_no_t,nDofsPerTargetElement> dofNosLocal;
      for (int dofIndex = 0; dofIndex < nDofsPerTargetElement; dofIndex++)
      {
        dofNosLocal[dofIndex] = fieldVariableTarget.functionSpace()->getDofNo(targetElementNoLocal, dofIndex);
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

#ifndef NDEBUG
  VecD<nComponents> previousTargetValue;
  VecD<FunctionSpaceSourceType::dim()> previousXi;
  element_no_t previousSourceElementNoLocal;
  std::array<double,FunctionSpaceSourceType::nDofsPerElement()> previousScalingFactors;
#endif

  // visualization for 1D-1D: s=source, t=target
  // s--t--------s-----t-----s

  // loop over all local dofs of the target functionSpace, which was the source function space when the mapping was initialized
  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    // if source dof is outside of target mesh
    if (!targetMappingInfo_[targetDofNoLocal].mapThisDof)
      continue;

    // the first set of surrounding nodes (targetElements[0]) is enough
    const typename targetDof_t::element_t &targetElement = targetMappingInfo_[targetDofNoLocal].targetElements[0];

    element_no_t sourceElementNoLocal = targetElement.elementNoLocal;

    // get source values of the element where targetDofNoLocal is in
    std::array<VecD<nComponents>,nDofsPerSourceElement> sourceValues;
    fieldVariableSource.getElementValues(sourceElementNoLocal, sourceValues);

    double scalingFactorsSum = 0;
    for (int i = 0; i < nDofsPerSourceElement; i++)
    {
      scalingFactorsSum += targetElement.scalingFactors[i];
    }

    VecD<nComponents> targetValue = sourceValues * targetElement.scalingFactors / scalingFactorsSum;
    fieldVariableTarget.setValue(targetDofNoLocal, targetValue);

#ifndef NDEBUG
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

    // the first set of surrounding nodes (targetElements[0]) is enough
    const typename targetDof_t::element_t &targetElement = targetMappingInfo_[targetDofNoLocal].targetElements[0];

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
