#include "mesh/mapping_between_meshes/manager/01_manager_initialize.h"

#include "function_space/function_space.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/unstructured_deformable.h"
#include "utility/python_utility.h"

namespace MappingBetweenMeshes
{

ManagerInitialize::ManagerInitialize(PythonConfig specificSettings) :
  ManagerLog(specificSettings)
{
  // parse all settings in "MappingsBetweenMeshes" and store them in mappingsBetweenMeshes_
  storeMappingsBetweenMeshes(specificSettings);
}

void ManagerInitialize::storeMappingBetweenMeshes(std::string sourceMeshName, PyObject *targetMeshPy)
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
      mappingWithSettings.xiTolerance = 0.1;
      mappingWithSettings.enableWarnings = false;
      mappingWithSettings.compositeUseOnlyInitializedMappings = false;
      mappingWithSettings.isEnabledFixUnmappedDofs = true;
      mappingWithSettings.defaultValue = 0.0;
      mappingsBetweenMeshes_[sourceMeshName].insert(std::pair<std::string,MappingWithSettings>(targetMeshToMapTo,mappingWithSettings));

      // log event, to be included in the log file
      addLogEntryMapping(sourceMeshName, targetMeshToMapTo, mappingLogEntry_t::logEvent_t::eventParseSettings);
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
    bool isEnabledFixUnmappedDofs = PythonUtility::getOptionBool(targetMeshPy, "fixUnmappedDofs", stringPath.str(), true);
    double defaultValue = PythonUtility::getOptionDouble(targetMeshPy, "defaultValue", stringPath.str(), 0.0);

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
      mappingWithSettings.isEnabledFixUnmappedDofs = isEnabledFixUnmappedDofs;
      mappingWithSettings.defaultValue = defaultValue;
      mappingsBetweenMeshes_[sourceMeshName].insert(std::pair<std::string,MappingWithSettings>(targetMeshToMapTo,mappingWithSettings));

      // log event, to be included in the log file
      addLogEntryMapping(sourceMeshName, targetMeshToMapTo, mappingLogEntry_t::logEvent_t::eventParseSettings);
    }
  }
  else
  {
    LOG(WARNING) << "Value for MappingsBetweenMeshes from mesh \"" << sourceMeshName << "\" should be either a string (the name of the mesh to map to)"
      << " or a dict, e.g. {\"name\": ..., \"xiTolerance\": 0.1, \"enableWarnings\": True, \"compositeUseOnlyInitializedMappings\": False} ";
  }
}

void ManagerInitialize::storeMappingsBetweenMeshes(PythonConfig specificSettings)
{
  LOG(TRACE) << "MeshManagerInitialize::storeMappingsBetweenMeshes";

  // retrieve python settings object
  if (specificSettings.pyObject())
  {
    // if the python setting contain the key "MappingsBetweenMeshes"
    std::string keyString("MappingsBetweenMeshes");
    if (specificSettings.hasKey(keyString))
    {
      // loop over the dict
      std::pair<std::string, PyObject *> dictItem
        = specificSettings.getOptionDictBegin<std::string, PyObject *>(keyString);

      for (; !specificSettings.getOptionDictEnd(keyString);
          specificSettings.getOptionDictNext<std::string, PyObject *>(keyString, dictItem))
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
  else 
  {
    LOG(WARNING) << "specificSettings is not set, this is probably in a unit test";
  }
}

} // namespace
