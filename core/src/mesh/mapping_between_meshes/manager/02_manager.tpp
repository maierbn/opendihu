#include "mesh/mapping_between_meshes/manager/02_manager.h"

namespace MappingBetweenMeshes
{

template<typename FunctionSpace1Type, typename FunctionSpace2Type>
void Manager::initializeMappingsBetweenMeshes(const std::shared_ptr<FunctionSpace1Type> functionSpace1, const std::shared_ptr<FunctionSpace2Type> functionSpace2)
{
  LOG(DEBUG) << "initializeMappingsBetweenMeshes source: \"" << functionSpace1->meshName() << "\" target: \"" << functionSpace2->meshName() << "\".";

  assert(FunctionSpace1Type::dim() <= FunctionSpace2Type::dim());   // the first mesh has to be lower dimensional than the second (or equal).

  // check if mapping functionSpace1 -> functionSpace2 is defined in config
  std::string sourceMeshName = functionSpace1->meshName();
  std::string targetMeshName = functionSpace2->meshName();

  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    LOG(DEBUG) << "key \"" << functionSpace1->meshName() << "\" exists";
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      if (mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping == nullptr)
      {
        LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" -> \"" << functionSpace2->meshName() << "\" declared, but not yet created, create mapping";

        // create the mapping from functionSpace1 to functionSpace2
        createMappingBetweenMeshes(functionSpace1, functionSpace2);
      }
      else
      {
        LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" -> \"" << functionSpace2->meshName() << "\" exists.";
      }
    }
    else
    {
      LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" -> \"" << functionSpace2->meshName() << "\" not declared";
    }
  }
  else
  {
    LOG(DEBUG) << "key \"" << functionSpace1->meshName() << "\" does not exist";
  }

  // change source and target mesh and initialize the reverse direction
  sourceMeshName = functionSpace2->meshName();
  targetMeshName = functionSpace1->meshName();

  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    LOG(DEBUG) << "key \"" << functionSpace2->meshName() << "\" exists";
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      if (mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping == nullptr)
      {
        LOG(DEBUG) << "\"" << functionSpace2->meshName() << "\" -> \"" << functionSpace1->meshName() << "\" declared, but not yet created, create mapping";

        // create the mapping from functionSpace2 to functionSpace1
        createMappingBetweenMeshes(functionSpace2, functionSpace1);
      }
      else
      {
        LOG(DEBUG) << "\"" << functionSpace1->meshName() << "\" -> \"" << functionSpace2->meshName() << "\" exists.";
      }
    }
    else
    {
      LOG(DEBUG) << "\"" << functionSpace2->meshName() << "\" -> \"" << functionSpace1->meshName() << "\" not declared";
    }
  }
  else
  {
    LOG(DEBUG) << "key \"" << functionSpace2->meshName() << "\" does not exist";
  }
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
bool Manager::
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

}   // namespace
