#include "mesh/mapping_between_meshes/manager/01_manager_initialize.h"

namespace MappingBetweenMeshes
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
bool ManagerInitialize::
hasMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                        std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget)
{
  assert(functionSpaceSource);
  assert(functionSpaceTarget);

  // get mesh names
  std::string sourceMeshName = functionSpaceSource->meshName();
  std::string targetMeshName = functionSpaceTarget->meshName();

  // check if the mapping exist
  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      return true;
    }
  }
  return false;
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>> ManagerInitialize::
createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget)
{
  std::string sourceMeshName = functionSpaceSource->meshName();
  std::string targetMeshName = functionSpaceTarget->meshName();

  if (!hasMappingBetweenMeshes(functionSpaceSource, functionSpaceTarget))
  {
    this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].xiTolerance = 0.1;
    this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].enableWarnings = false;
    this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].compositeUseOnlyInitializedMappings = false;
    this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].isEnabledFixUnmappedDofs = false;
  }
  else if (this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping)
  {
    // check if the mapping already existed
    LOG(WARNING) << "Mapping from mesh \"" << sourceMeshName << "\" to mesh \"" << targetMeshName << "\" is already defined.";
  }

  // get options for the mapping
  double xiTolerance = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].xiTolerance;
  bool enableWarnings = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].enableWarnings;
  bool compositeUseOnlyInitializedMappings = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].compositeUseOnlyInitializedMappings;
  bool isEnabledFixUnmappedDofs = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].isEnabledFixUnmappedDofs;

  LOG(INFO) << "create MappingBetweenMeshes \"" << sourceMeshName << "\" (" << FunctionSpaceSourceType::dim() << "D, " << functionSpaceSource->nNodesGlobal() << " nodes) -> \""
     << targetMeshName << "\" (" << FunctionSpaceTargetType::dim() << "D, " << functionSpaceTarget->nNodesGlobal() << " nodes), xiTolerance: " << xiTolerance;

  // log event, to be included in the log file
  addLogEntryMapping(functionSpaceSource, functionSpaceTarget, mappingLogEntry_t::logEvent_t::eventCreateMapping);

  // create the mapping under the given source and target mesh names
  this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping = std::static_pointer_cast<MappingBetweenMeshesBase>(
    std::make_shared<MappingBetweenMeshes<FunctionSpaceSourceType,FunctionSpaceTargetType>>(functionSpaceSource, functionSpaceTarget,
                                                                                            xiTolerance, enableWarnings, compositeUseOnlyInitializedMappings, isEnabledFixUnmappedDofs)
  );

  return std::static_pointer_cast<MappingBetweenMeshes<FunctionSpaceSourceType,FunctionSpaceTargetType>>(
    this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping);
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
std::shared_ptr<MappingBetweenMeshes<typename FunctionSpaceSourceType::FunctionSpace, typename FunctionSpaceTargetType::FunctionSpace>> ManagerInitialize::
mappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                     std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget)
{
  typedef MappingBetweenMeshes<typename FunctionSpaceSourceType::FunctionSpace, typename FunctionSpaceTargetType::FunctionSpace> MappingType;

  assert(functionSpaceSource);
  assert(functionSpaceTarget);

  // get mesh names
  std::string sourceMeshName = functionSpaceSource->meshName();
  std::string targetMeshName = functionSpaceTarget->meshName();

  // check if the mapping exists already
  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      std::shared_ptr<MappingBetweenMeshesBase> mappingBase = mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping;

      // if mapping has already been created, return it
      if (mappingBase)
      {
        return std::static_pointer_cast<MappingType>(mappingBase);
      }
    }
  }

  // mapping between meshes has not yet been created, create now
  // if it does not yet exist, output message and create it
  LOG(DEBUG) << "Mapping from mesh \"" << sourceMeshName << "\" to \"" << targetMeshName
    << "\" was not initialized. Initializing now. Specify MappingsBetweenMeshes { \"" << sourceMeshName << "\" : \"" << targetMeshName << "\" } as top level object of the python config. "
    << "(It could be that this was done, but because of MultipleInstances in the OperatorSplitting or Coupling, only the first mesh mapping got initialized.)";

  // create the mapping
  std::shared_ptr<MappingType> mapping;
  mapping = createMappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>(
    functionSpaceSource, functionSpaceTarget
  );
  return mapping;
}

template<typename FunctionSpace1Type, typename FunctionSpace2Type>
void ManagerInitialize::
initializeMappingsBetweenMeshesFromSettings(const std::shared_ptr<FunctionSpace1Type> functionSpace1,
                                            const std::shared_ptr<FunctionSpace2Type> functionSpace2)
{
  // initialize functionSpace1 <-> functionSpace2 in both directions, only the mapping that have been defined in the settings will be created
  if (hasMappingBetweenMeshes(functionSpace1, functionSpace2))
  {
    mappingBetweenMeshes<FunctionSpace1Type,FunctionSpace2Type>(functionSpace1,functionSpace2);

    if (hasMappingBetweenMeshes(functionSpace2, functionSpace1))
    {
      LOG(WARNING) << "A mapping between the meshes is defined in both directions. "
        << "There exist settings for \"" << functionSpace1->meshName() << "\" -> " << functionSpace2->meshName()
        << "\" as well as for \"" << functionSpace2->meshName() << "\" -> " << functionSpace1->meshName() << "\".\n"
        << "Usually it makes sense to only define the mapping:\n   - from the lower-dimensional mesh to the higher-dimensional mesh\n   "
        << "- or for meshes with equal dimensions from the fine mesh to the coarse mesh.";

      mappingBetweenMeshes<FunctionSpace2Type,FunctionSpace1Type>(functionSpace2,functionSpace1);
    }
  }
  else
  {
    if (hasMappingBetweenMeshes(functionSpace2, functionSpace1))
    {
      mappingBetweenMeshes<FunctionSpace2Type,FunctionSpace1Type>(functionSpace2,functionSpace1);
    }
    else
    {
      // there is no mapping defined between functionSpace1 and functionSpace2
    }
  }
}


}   // namespace
