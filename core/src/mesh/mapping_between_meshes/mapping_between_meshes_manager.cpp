#include "mesh/mapping_between_meshes/mapping_between_meshes_manager.h"

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
        else if (PyUnicode_Check(value))
        {
          std::string targetMeshToMapTo = PythonUtility::convertFromPython<std::string>::get(value);
          VLOG(1) << "Store mapping between mesh \"" << key << "\" and " << targetMeshToMapTo;

          if (mappingsBetweenMeshes_[key].find(key) == mappingsBetweenMeshes_[key].end())
          {
            MappingWithSettings mappingWithSettings;
            mappingWithSettings.mapping = nullptr;
            mappingWithSettings.xiTolerance = 0.0;
            mappingsBetweenMeshes_[key].insert(std::pair<std::string,MappingWithSettings>(targetMeshToMapTo,mappingWithSettings));
          }
        }
        else if (PyDict_Check(value))
        {
          std::stringstream stringPath;
          stringPath << specificSettings_.getStringPath();
          stringPath << "[\"" << value << "\"]";
          std::string targetMeshToMapTo = PythonUtility::getOptionString(value, "name", stringPath.str(), "");
          double xiTolerance = PythonUtility::getOptionDouble(value, "xiTolerance", stringPath.str(), 0.1);    // 0.1 was tested to be reasonable, the actual tolerance in computation of the xi values of points inside elements is independent of this value

          VLOG(1) << "Store mapping between mesh \"" << key << "\" and " << targetMeshToMapTo << " with xiTolerance " << xiTolerance;

          if (mappingsBetweenMeshes_[key].find(key) == mappingsBetweenMeshes_[key].end())
          {
            MappingWithSettings mappingWithSettings;
            mappingWithSettings.mapping = nullptr;
            mappingWithSettings.xiTolerance = xiTolerance;
            mappingsBetweenMeshes_[key].insert(std::pair<std::string,MappingWithSettings>(targetMeshToMapTo,mappingWithSettings));
          }
        }
        else
        {
          LOG(WARNING) << "Value for MappingsBetweenMeshes from mesh \"" << key << "\" should be either a string (the name of the mesh to map to)"
            << " or a dict {\"name\": ..., \"xiTolerance\": ...} ";
        }
      }
    }
  }
}

} // namespace
