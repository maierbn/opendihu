#include "mesh/mapping_between_meshes_manager.h"

#include "function_space/function_space.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/unstructured_deformable.h"

namespace Mesh
{

MappingBetweenMeshesManager::MappingBetweenMeshesManager(PythonConfig specificSettings) :
  specificSettings_(specificSettings)
{
  LOG(TRACE) << "MeshMappingBetweenMeshesManager constructor";
  storeMappingsBetweenMeshes();
}

void MappingBetweenMeshesManager::storeMappingsBetweenMeshes()
{
  LOG(TRACE) << "MeshMappingBetweenMeshesManager::storeMappingsBetweenMeshes";
  if (specificSettings_.pyObject())
  {
    std::string keyString("MappingsBetweenMeshes");
    if (specificSettings_.hasKey(keyString))
    {
      std::pair<std::string, PyObject *> dictItem
        = specificSettings_.getOptionDictBegin<std::string, PyObject *>(keyString);

      for (; !specificSettings_.getOptionDictEnd(keyString);
          specificSettings_.getOptionDictNext<std::string, PyObject *>(keyString, dictItem))
      {
        std::string key = dictItem.first;
        PyObject *value = dictItem.second;

        if (value == NULL)
        {
          LOG(WARNING) << "Could not extract dict for MappingsBetweenMeshes[\"" << key << "\"].";
        }
        else if (!PyUnicode_Check(value))
        {
          LOG(WARNING) << "Value for MappingsBetweenMeshes from mesh \"" << key << "\" should be a string (the name of the mesh to map to).";
        }
        else
        {
          std::string targetMeshToMapTo = PythonUtility::convertFromPython<std::string>::get(value);
          LOG(DEBUG) << "Store mapping between mesh \"" << key << "\" and " << targetMeshToMapTo;

          if (mappingsBetweenMeshes_[key].find(key) == mappingsBetweenMeshes_[key].end())
          {
            mappingsBetweenMeshes_[key].insert(std::pair<std::string,std::shared_ptr<MappingBetweenMeshesBase>>(targetMeshToMapTo,nullptr));
          }
        }
      }
    }
  }
}

std::shared_ptr<MappingBetweenMeshesBase> MappingBetweenMeshesManager::mappingBetweenMeshes(std::string sourceMeshName, std::string targetMeshName)
{
  if (mappingsBetweenMeshes_.find(sourceMeshName) == mappingsBetweenMeshes_.end())
    return nullptr;

  if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) == mappingsBetweenMeshes_[sourceMeshName].end())
    return nullptr;

  return mappingsBetweenMeshes_[sourceMeshName][targetMeshName];
}

} // namespace
