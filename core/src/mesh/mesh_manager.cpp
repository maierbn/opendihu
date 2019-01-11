#include "mesh/mesh_manager.h"

#include "function_space/function_space.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/unstructured_deformable.h"

namespace Mesh
{

Manager::Manager(PythonConfig specificSettings) :
  partitionManager_(nullptr), specificSettings_(specificSettings), numberAnonymousMeshes_(0)
{
  LOG(TRACE) << "MeshManager constructor";
  storePreconfiguredMeshes();
  storeMappingsBetweenMeshes();
}

void Manager::setPartitionManager(std::shared_ptr<Partition::Manager> partitionManager)
{
  partitionManager_ = partitionManager;
}

void Manager::storePreconfiguredMeshes()
{
  LOG(TRACE) << "MeshManager::storePreconfiguredMeshes";
  if (specificSettings_.pyObject())
  {
    std::string keyString("Meshes");
    if (specificSettings_.hasKey("Meshes"))
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
          LOG(WARNING) << "Could not extract dict for Mesh \"" << key << "\".";
        }
        else if (!PyDict_Check(value))
        {
          LOG(WARNING) << "Value for mesh with name \"" << key << "\" should be a dict.";
        }
        else
        {
          LOG(DEBUG) << "Store mesh configuration with key \"" << key << "\".";
          if (meshConfiguration_.find(key) != meshConfiguration_.end())
          {
            meshConfiguration_.at(key).setPyObject(value);
          }
          else
          {
            meshConfiguration_.insert(std::pair<std::string,PythonConfig>(key, PythonConfig(specificSettings_, "Meshes", key, value)));
          }
        }
      }
    }
    else
    {
      LOG(INFO) << "You have specified the mesh in-line and not under the extra key \"Meshes\". You could do so,"
        " by defining \"Meshes\": {\"<your custom mesh name>\": {<your mesh parameters>}} at the beginning of the"
        " config and \"meshName\": \"<your custom mesh name>\" where you currently have specified the mesh parameters."
        " This is required if you want to use the same mesh for multiple objects.";
    }
  }
}

void Manager::storeMappingsBetweenMeshes()
{
  LOG(TRACE) << "MeshManager::storeMappingsBetweenMeshes";
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

std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::Generic,1>> Manager::
createGenericFieldVariable(int nEntries, std::string name)
{
  assert(nEntries > 1);

  // create generic field variable

  // constructor is declared in function_space/06_function_space_dofs_nodes.h
  // FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent);

  std::array<element_no_t, 1> nElements({nEntries - 1});
  std::array<double, 1> physicalExtent({0.0});
  std::stringstream meshName;
  meshName << "meshForFieldVariable" << name;
  std::shared_ptr<Mesh> mesh = createFunctionSpace<FunctionSpace::Generic>(meshName.str(), nElements, physicalExtent);

  LOG(DEBUG) << "create generic field variable with " << nEntries << " entries.";
  std::shared_ptr<FunctionSpace::Generic> functionSpace = std::static_pointer_cast<FunctionSpace::Generic>(mesh);

  // createFieldVariable is declared in function_space/10_function_space_field_variable.h
  //template <int nComponents>
  //std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name);
  return functionSpace->template createFieldVariable<1>(name);
}

bool Manager::hasFunctionSpace(std::string meshName)
{
  return functionSpaces_.find(meshName) != functionSpaces_.end();
}

std::shared_ptr<MappingBetweenMeshesBase> Manager::mappingBetweenMeshes(std::string sourceMeshName, std::string targetMeshName)
{
  if (mappingsBetweenMeshes_.find(sourceMeshName) == mappingsBetweenMeshes_.end())
    return nullptr;

  if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) == mappingsBetweenMeshes_[sourceMeshName].end())
    return nullptr;

  return mappingsBetweenMeshes_[sourceMeshName][targetMeshName];
}

} // namespace
