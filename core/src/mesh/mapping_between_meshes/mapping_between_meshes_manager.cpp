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

void MappingBetweenMeshesManager::storeMappingBetweenMeshes(std::string sourceMeshName, PyObject *targetMeshPy)
{
  if (PyUnicode_Check(targetMeshPy))
  {
    std::string targetMeshToMapTo = PythonUtility::convertFromPython<std::string>::get(targetMeshPy);
    VLOG(1) << "Store mapping between mesh \"" << sourceMeshName << "\" and " << targetMeshToMapTo;

    // if the mapping settings have not yet been stored
    if (mappingsBetweenMeshes_[sourceMeshName].find(sourceMeshName) == mappingsBetweenMeshes_[sourceMeshName].end())
    {
      // create new entry of settings in mappingsBetweenMeshes_
      MappingWithSettings mappingWithSettings;
      mappingWithSettings.mapping = nullptr;
      mappingWithSettings.xiTolerance = 0.0;
      mappingWithSettings.enableWarnings = true;
      mappingWithSettings.compositeUseOnlyInitializedMappings = false;
      mappingsBetweenMeshes_[sourceMeshName].insert(std::pair<std::string,MappingWithSettings>(targetMeshToMapTo,mappingWithSettings));
    }
  }
  else if (PyDict_Check(targetMeshPy))
  {
    std::stringstream stringPath;
    stringPath << "[\"MappingsBetweenMeshes\"][\"" << sourceMeshName << "\"]";

    // parse the name of the target mesh
    std::string targetMeshToMapTo = PythonUtility::getOptionString(targetMeshPy, "name", stringPath.str(), "");

    // parse the xi tolerance of the mapping
    double xiTolerance = PythonUtility::getOptionDouble(targetMeshPy, "xiTolerance", stringPath.str(), 0.1);    // 0.1 was tested to be reasonable, the actual tolerance in computation of the xi values of points inside elements is independent of this value

    // parse if warnings should be shown
    bool enableWarnings = PythonUtility::getOptionBool(targetMeshPy, "enableWarnings", stringPath.str(), true);
    bool compositeUseOnlyInitializedMappings = PythonUtility::getOptionBool(targetMeshPy, "compositeUseOnlyInitializedMappings", stringPath.str(), false);

    VLOG(1) << "Store mapping between mesh \"" << sourceMeshName << "\" and " << targetMeshToMapTo << " with xiTolerance " << xiTolerance;

    // if the mapping settings have not yet been stored
    if (mappingsBetweenMeshes_[sourceMeshName].find(sourceMeshName) == mappingsBetweenMeshes_[sourceMeshName].end())
    {
      // store the new settings
      MappingWithSettings mappingWithSettings;
      mappingWithSettings.mapping = nullptr;
      mappingWithSettings.xiTolerance = xiTolerance;
      mappingWithSettings.enableWarnings = enableWarnings;
      mappingWithSettings.compositeUseOnlyInitializedMappings = compositeUseOnlyInitializedMappings;
      mappingsBetweenMeshes_[sourceMeshName].insert(std::pair<std::string,MappingWithSettings>(targetMeshToMapTo,mappingWithSettings));
    }
  }
  else
  {
    LOG(WARNING) << "Value for MappingsBetweenMeshes from mesh \"" << sourceMeshName << "\" should be either a string (the name of the mesh to map to)"
      << " or a dict, e.g. {\"name\": ..., \"xiTolerance\": 0.1, \"enableWarnings\": True, \"compositeUseOnlyInitializedMappings\": False} ";
  }
}

void MappingBetweenMeshesManager::storeMappingsBetweenMeshes()
{
  LOG(TRACE) << "MeshMappingBetweenMeshesManager::storeMappingsBetweenMeshes";

  // retrieve python settings object
  if (specificSettings_.pyObject())
  {
    // if the python setting contain the key "MappingsBetweenMeshes"
    std::string keyString("MappingsBetweenMeshes");
    if (specificSettings_.hasKey(keyString))
    {
      // loop over the dict
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
        else if (PyList_Check(value))
        {
          // value is a list of multiple meshes
          std::vector<PyObject *> targetMeshesPy = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(value);
          for (PyObject *targetMeshPy : targetMeshesPy)
          {
            storeMappingBetweenMeshes(key, targetMeshPy);
          }
        }
        else
        {
          storeMappingBetweenMeshes(key, value);
        }
      }
    }
  }
}

} // namespace
