#include "mesh/mapping_between_meshes.h"

namespace Mesh
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                                                                                             std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget) :
  functionSpaceSource_(functionSpaceSource),
  functionSpaceTarget_(functionSpaceTarget)
{

  const dof_no_t nDofsLocalSource = functionSpaceSource->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  element_no_t elementNo;
  int ghostMeshNo;
  std::array<double,FunctionSpaceTargetType::dim()> xi;
  bool startSearchInCurrentElement = false;

  targetMappingInfo_.resize(nDofsLocalSource);

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
    // determine information how to map a source value to the target mesh
    targetDof_t targetMappingInfo;

    // get node position of the source dof
    Vec3 position;
    dof_no_t sourceDofNoGlobal = functionSpaceTarget->meshPartition()->getDofNoGlobalPetsc(sourceDofNoLocal);
    functionSpaceTarget->getGeometry(sourceDofNoGlobal);

    // find element no in the target mesh where the position is
    if (!functionSpaceTarget->findPosition(position, elementNo, ghostMeshNo, xi, startSearchInCurrentElement))
    {
      LOG(ERROR) << "Could not create mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
        << functionSpaceTarget->meshName() << "\", dof local " << sourceDofNoLocal << ", global " << sourceDofNoGlobal
        << " of mesh \"" << functionSpaceSource->meshName() << "\" at position " << position << " is outside of mesh \""
        << functionSpaceTarget->meshName() << "\".";
        break;
    }

    // store element no
    targetMappingInfo.elementNoLocal = elementNo;

    // determine factors how to distribute the value to the dofs of the target element
    targetMappingInfo.scalingFactors.resize(nDofsPerTargetElement);

    // note: geometry value = sum over dofs of geometryValue_dof * phi_dof(xi)
    for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
    {
      targetMappingInfo.scalingFactors[targetDofIndex] = functionSpaceTarget->phi(targetDofIndex,xi);
    }

    targetMappingInfo_[sourceDofNoLocal] = targetMappingInfo;

    // next time when searching for the target element, start search from previous element
    startSearchInCurrentElement = true;
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  template<int nComponentsSource, int nComponentsTarget>
void MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>::map(
  FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableSource, int componentNoSource,
  FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableTarget, int componentNoTarget)
{
  assert(componentNoSource >= 0 && componentNoSource < nComponentsSource);
  assert(componentNoTarget >= 0 && componentNoTarget < nComponentsTarget);

  fieldVariableTarget->zeroValues();

  const dof_no_t nDofsLocalSource = fieldVariableSource->functionSpace()->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  // loop over all local dofs of the source functionSpace
  for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
  {
    // get value
    double sourceValue;
    fieldVariableSource->getValue(componentNoSource, sourceDofNoLocal, sourceValue);

    // store the value to the target function space
    element_no_t targetElementNoLocal = targetMappingInfo_[sourceDofNoLocal].elementNoLocal;
    std::array<double,nDofsPerTargetElement> targetValues = targetMappingInfo_[sourceDofNoLocal].scalingFactors * sourceValue;
    fieldVariableTarget->setElementalValues(componentNoTarget, targetElementNoLocal, targetValues, ADD_VALUES);
  }
}

}  // namespace
