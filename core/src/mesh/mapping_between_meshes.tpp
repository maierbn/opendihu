#include "mesh/mapping_between_meshes.h"

#include "control/performance_measurement.h"

namespace Mesh
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                                                                                             std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget) :
  functionSpaceSource_(functionSpaceSource),
  functionSpaceTarget_(functionSpaceTarget)
{
  // create the mapping

  Control::PerformanceMeasurement::start("compute mapping");

  const dof_no_t nDofsLocalSource = functionSpaceSource->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  element_no_t elementNo = 0;
  int ghostMeshNo = 0;
  std::array<double,FunctionSpaceTargetType::dim()> xi;

  targetMappingInfo_.resize(nDofsLocalSource);

  VLOG(1) << "create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();

  //VLOG(1) << "target meshPartition: " << *functionSpaceTarget->meshPartition();
  //VLOG(1) << "geometryField: " << functionSpaceTarget->geometryField();

  double xiToleranceBase = 1e-2;

  bool mappingSucceeded = true;
  bool startSearchInCurrentElement = false;

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
    for(nTries = 0; nTries < 10; nTries++)
    {
      // find element no in the target mesh where the position is
      if (functionSpaceTarget->findPosition(position, elementNo, ghostMeshNo, xi, startSearchInCurrentElement, xiTolerance))
      {
        xiToleranceBase = xiTolerance;
        break;
      }
      else
      {
        xiTolerance *= 2;
        LOG(DEBUG) << "Try again with xiTolerance = " << xiTolerance;
      }
    }

    if (nTries == 10)
    {
      LOG(DEBUG) << "Could not create mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
        << functionSpaceTarget->meshName() << "\", dof local " << sourceDofNoLocal
        << " of mesh \"" << functionSpaceSource->meshName() << "\" at position " << position << " is outside of mesh \""
        << functionSpaceTarget->meshName() << "\" with tolerance " << xiTolerance << ".";
      mappingSucceeded = false;
      break;
    }

    // store element no
    targetMappingInfo.elementNoLocal = elementNo;

    // determine factors how to distribute the value to the dofs of the target element

    // note: geometry value = sum over dofs of geometryValue_dof * phi_dof(xi)
    for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
    {
      targetMappingInfo.scalingFactors[targetDofIndex] = functionSpaceTarget->phi(targetDofIndex, xi);
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

  Control::PerformanceMeasurement::stop("compute mapping");

  if (!mappingSucceeded)
  {
    LOG(ERROR) << "Could not create mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
      << functionSpaceTarget->meshName() << "\".";
    LOG(FATAL) << "end";
  }
  else
  {
    LOG(DEBUG) << "Successfully initialized mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
      << functionSpaceTarget->meshName() << "\".";
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  template<int nComponentsSource, int nComponentsTarget>
void MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::map(
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

  //std::vector<dof_no_t> sourceLocalDofNos(nDofsLocalSource);
  //std::iota(sourceLocalDofNos.begin(), sourceLocalDofNos.end(), 0);
  //fieldVariableSource.getValues(componentNoSource, sourceLocalDofNos, sourceValues);

  VLOG(1) << "map " << fieldVariableSource.name() << "." << componentNoSource <<
    " (" << fieldVariableSource.functionSpace()->meshName() << ") -> " << fieldVariableTarget.name() << "." << componentNoTarget
    << " (" << fieldVariableTarget.functionSpace()->meshName() << ")";

  VLOG(1) << "source has " << nDofsLocalSource << " local dofs";
  VLOG(1) << fieldVariableSource;
  VLOG(1) << "extracted source values: " << sourceValues;

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
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

    VLOG(2) << "  source dof " << sourceDofNoLocal << ", value: " << sourceValue << ", scaling factors: " << targetMappingInfo_[sourceDofNoLocal].scalingFactors
      << ", targetValues: " << targetValues
      << ", targetElementNoLocal: " << targetElementNoLocal << ", target dofs: " << dofNosLocal;
  }
}

}  // namespace
