#include "mesh/mapping_between_meshes_manager.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/structured_regular_fixed.h"
#include "control/performance_measurement.h"

namespace Mesh
{

template<typename FunctionSpace1Type, typename FunctionSpace2Type>
void MappingBetweenMeshesManager::initializeMappingsBetweenMeshes(const std::shared_ptr<FunctionSpace1Type> functionSpace1, const std::shared_ptr<FunctionSpace2Type> functionSpace2)
{
  LOG(DEBUG) << "initializeMappingsBetweenMeshes source: \"" << functionSpace1->meshName() << "\" target: \"" << functionSpace2->meshName() << "\".";

  // check if mapping functionSpace1 -> functionSpace2 is defined in config
  //assert(functionSpaces_.find(functionSpace1->meshName()) != functionSpaces_.end());
  //assert(functionSpaces_.find(functionSpace2->meshName()) != functionSpaces_.end());

  std::string sourceMeshName = functionSpace1->meshName();
  std::string targetMeshName = functionSpace2->meshName();

  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" declared";
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" -> \"" << functionSpace2->meshName() << "\" declared";
      mappingsBetweenMeshes_[sourceMeshName][targetMeshName] = std::static_pointer_cast<MappingBetweenMeshesBase>(
        std::make_shared<MappingBetweenMeshes<FunctionSpace1Type,FunctionSpace2Type>>(functionSpace1,functionSpace2)
      );
    }
    else
    {
      LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" -> \"" << functionSpace2->meshName() << "\" not declared";
    }
  }
  else
  {
    LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" not declared";
  }

  sourceMeshName = functionSpace2->meshName();
  targetMeshName = functionSpace1->meshName();

  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    LOG(DEBUG) << "\"" << functionSpace2->meshName() << "\" declared";
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      LOG(DEBUG) << "\"" << functionSpace2->meshName() << "\" -> \"" << functionSpace1->meshName() << "\" declared";
      mappingsBetweenMeshes_[sourceMeshName][targetMeshName] = std::static_pointer_cast<MappingBetweenMeshesBase>(
        std::make_shared<MappingBetweenMeshes<FunctionSpace2Type,FunctionSpace1Type>>(functionSpace2,functionSpace1)
      );
    }
    else
    {
      LOG(DEBUG) << "\"" << functionSpace2->meshName() << "\" -> \"" << functionSpace1->meshName() << "\" not declared";
    }
  }
  else
  {
    LOG(DEBUG) << "\"" << functionSpace2->meshName() << "\" not declared";
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>> MappingBetweenMeshesManager::
createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget)
{
  std::string sourceMeshName = functionSpaceSource->meshName();
  std::string targetMeshName = functionSpaceTarget->meshName();

  if (this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName])
  {
    LOG(WARNING) << "Mapping from mesh \"" << sourceMeshName << "\" to mesh \"" << targetMeshName << "\" is already defined.";
  }

  this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName] = std::static_pointer_cast<MappingBetweenMeshesBase>(
    std::make_shared<MappingBetweenMeshes<FunctionSpaceSourceType,FunctionSpaceTargetType>>(functionSpaceSource,functionSpaceTarget)
  );

  return std::static_pointer_cast<MappingBetweenMeshes<FunctionSpaceSourceType,FunctionSpaceTargetType>>(this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName]);
}

//! prepare a mapping to the fieldVariableTarget, this zeros the targetFactorSum for the field variable
template<typename FieldVariableTargetType>
void MappingBetweenMeshesManager::
prepareMapping(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
{
  std::string fieldVariableTargetName = fieldVariableTarget->name();
  // std::map<std::string, std::shared_ptr<FieldVariable::FieldVariableBase>> targetFactorSum_;

  // if the targetFactorSum field variable does not yet exist, create it
  if (targetFactorSum_.find(fieldVariableTargetName) == targetFactorSum_.end())
  {
    std::stringstream name;
    name << "targetFactorSum_" << fieldVariableTargetName;

    std::vector<std::string> componentNames(1, "0");

    targetFactorSum_[fieldVariableTargetName] = std::static_pointer_cast<FieldVariable::FieldVariableBase>(
      std::make_shared<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(*fieldVariableTarget, name.str(), componentNames)
    );
  }

  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(targetFactorSum_[fieldVariableTarget->name()]);

  // set all entries to 0
  targetFactorSum->zeroEntries();
  targetFactorSum->zeroGhostBuffer();
  fieldVariableTarget->zeroEntries();
  fieldVariableTarget->zeroGhostBuffer();
}

//! map data from the source to the target field variable. This has to be called between prepareMapping and finalizeMapping, can be called multiple times with different source meshes.
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void MappingBetweenMeshesManager::
map(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
    std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  std::string sourceMeshName = fieldVariableSource->functionSpace()->meshName();
  std::string targetMeshName = fieldVariableTarget->functionSpace()->meshName();

  typedef MappingBetweenMeshes<typename FieldVariableSourceType::FunctionSpace, typename FieldVariableTargetType::FunctionSpace> MappingType;

  std::shared_ptr<MappingType> mapping = nullptr;

  // try to get the stored mapping if it exists
  std::shared_ptr<MappingBetweenMeshesBase> mappingBase = mappingBetweenMeshes(sourceMeshName, targetMeshName);

  // if it exists, convert back to MappingType
  if (mappingBase)
  {
    mapping = std::static_pointer_cast<MappingType>(mappingBase);
  }
  else
  {
    // if it does not yet exist, output message and create it
    LOG(DEBUG) << "Mapping from mesh \"" << sourceMeshName << "\" to \"" << targetMeshName
      << "\" was not initialized. Initializing now. Specify MappingsBetweenMeshes { \"" << sourceMeshName << "\" : \"" << targetMeshName << "\" } as top level object of the python config. "
      << "(It could be that this was done, but because of MultipleInstances in the OperatorSplitting or Coupling, only the first mesh mapping got initialized.)";

    mapping = createMappingBetweenMeshes<typename FieldVariableSourceType::FunctionSpace, typename FieldVariableTargetType::FunctionSpace>(
      fieldVariableSource->functionSpace(), fieldVariableTarget->functionSpace()
    );
  }

  Control::PerformanceMeasurement::start("map");

  // assert that targetFactorSum_ field variable exists, this should have been created by prepareMapping()
  assert(targetFactorSum_.find(fieldVariableTarget->name()) != targetFactorSum_.end());

  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(targetFactorSum_[fieldVariableTarget->name()]);

  // call the method of the mapping that does not actual data transfer
  mapping->template map<FieldVariableSourceType::nComponents(),FieldVariableTargetType::nComponents()>(
    *fieldVariableSource, componentNoSource, *fieldVariableTarget, componentNoTarget, *targetFactorSum
  );

  Control::PerformanceMeasurement::stop("map");
}

//! finalize the mapping to the fieldVariableTarget, this computes the final values at the dofs from the accumulated values by dividing by the targetFactorSums
template<typename FieldVariableTargetType>
void MappingBetweenMeshesManager::
finalizeMapping(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  // assert that targetFactorSum_ field variable exists, this should have been created by prepareMapping()
  assert(targetFactorSum_.find(fieldVariableTarget->name()) != targetFactorSum_.end());


  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(targetFactorSum_[fieldVariableTarget->name()]);

  const dof_no_t nDofsLocalTarget = fieldVariableTarget->nDofsLocalWithoutGhosts();

  // communicate ghost values
  fieldVariableTarget->finishGhostManipulation();
  targetFactorSum->finishGhostManipulation();

  // compute final values by dividing by factorSum
  std::vector<double> targetValues;
  fieldVariableTarget->getValuesWithoutGhosts(componentNoTarget, targetValues);

  std::vector<double> targetFactorSums;
  targetFactorSum->getValuesWithoutGhosts(targetFactorSums);

  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    VLOG(2) << "  target dof " << targetDofNoLocal << ", divide value " << targetValues[targetDofNoLocal]
      << " by " << targetFactorSums[targetDofNoLocal] << ": " << targetValues[targetDofNoLocal]/targetFactorSums[targetDofNoLocal];
    if (fabs(targetFactorSums[targetDofNoLocal]) > 1e-12)
    {
      targetValues[targetDofNoLocal] /= targetFactorSums[targetDofNoLocal];
    }
    else
    {
      // output warning message, compute helper variables

      // get node no from dof no, this is needed for the global coordinates later
      node_no_t nodeNoLocal = int(targetDofNoLocal/fieldVariableTarget->functionSpace()->nDofsPerNode());

      // get current global node coordinates
      std::stringstream s;
      for (int i = 0; i < FieldVariableTargetType::FunctionSpace::dim(); i++)
      {
        if (i != 0)
        {
          s << ",";
        }
        s << fieldVariableTarget->functionSpace()->meshPartition()->nNodesGlobal(i);
      }

      // output the warning
      LOG(WARNING) << "In mapping to " << fieldVariableTarget->name() << "." << componentNoTarget
        << " (" << fieldVariableTarget->functionSpace()->meshName() << "), no values for target dof " << targetDofNoLocal
        << ", coordinates global: " << fieldVariableTarget->functionSpace()->meshPartition()->getCoordinatesGlobal(nodeNoLocal)
        << " of (" << s.str() << ")! Assuming 0.0.";
    }
  }

  // set the computed values
  fieldVariableTarget->setValuesWithoutGhosts(componentNoTarget, targetValues);
}


}   // namespace
