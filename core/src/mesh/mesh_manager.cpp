#include "mesh/mesh_manager.h"

MeshManager::MeshManager(const DihuContext& context) :
  context_(context), numberAnonymousMeshes_(0)
{
  LOG(TRACE) << "MeshManager constructor";
  specificSettings_ = this->context_.getPythonConfig();
  storePreconfiguredMeshes();
}

void MeshManager::storePreconfiguredMeshes()
{
  LOG(TRACE) << "MeshManager::storePreconfiguredMeshes";
  if (specificSettings_)
  {
    std::string keyString("Meshes");
    std::pair<std::string, PyObject *> dictItem 
      = PythonUtility::getOptionDictBegin<std::string, PyObject *>(specificSettings_, keyString);
    
    for (; !PythonUtility::getOptionDictEnd(specificSettings_, keyString); 
        PythonUtility::getOptionDictNext<std::string, PyObject *>(specificSettings_, keyString, dictItem))
    {
      std::string key = dictItem.first;
      PyObject *value = dictItem.second;
          
      if (value == NULL)
      {
        LOG(WARNING) << "Could not extract dict for Mesh \""<<key<<"\".";
      }
      else if(!PyDict_Check(value))
      {
        LOG(WARNING) << "Value for mesh with name \""<<key<<"\" should be a dict.";
      }
      else
      {
        LOG(DEBUG) << "store mesh configuration with key \""<<key<<"\".";
        meshConfiguration_[key] = value;
      }
    }
  }
}

bool MeshManager::hasMesh(std::string meshName)
{
  LOG(DEBUG) << "hasMesh("<<meshName<<")";
  LOG(DEBUG) << "meshes size: " << meshes_.size();
  
  return meshes_.find(meshName) != meshes_.end();
}

template<>
std::shared_ptr<Mesh::Mesh> MeshManager::mesh<Mesh::Mesh>(PyObject *settings)
{
  std::string meshName;
  if (PythonUtility::containsKey(settings, "meshName"))
  {
    meshName = PythonUtility::getOptionString(settings, "meshName", "");
    LOG(DEBUG) << "Config contains meshName \""<<meshName<<"\".";
    
    if (hasMesh(meshName))
    {
      return meshes_[meshName];
    }
    else if(meshConfiguration_.find(meshName) != meshConfiguration_.end())
    {
      // mesh was preconfigured, do nothing specific here, created standard mesh with 1 node
      LOG(DEBUG) << "Mesh configuration for \""<<meshName<<"\" found and requested, will be created now. "
        << " Type is not clear, so go for RegularFixed<1>.";
      meshes_[meshName] = std::make_shared<Mesh::RegularFixed<1>>(meshConfiguration_[meshName]);
      LOG(DEBUG) << "Stored under key \""<<meshName<<"\".";
      return std::static_pointer_cast<Mesh::Mesh>(meshes_[meshName]);
    }
    else
    {
      LOG(ERROR) << "Config contains reference to mesh with meshName \""<<meshName<<"\" but no such mesh was defined.";      
    }
  }
  else
  {
    LOG(DEBUG) << "Config does not contain meshName.";
  }
  
  // if nElements was specified, create standard regularFixed mesh
  if (PythonUtility::containsKey(settings, "nElements"))
  {
    // create new mesh, store as anonymous object
    std::stringstream anonymousName;
    anonymousName << "anonymous" << numberAnonymousMeshes_++;
    LOG(DEBUG) << "Create new mesh with type "<<typeid(Mesh::RegularFixed<1>).name()<<" and name \""<<anonymousName.str()<<"\".";
    std::shared_ptr<Mesh::Mesh> mesh = std::make_shared<Mesh::RegularFixed<1>>(settings);
    meshes_[anonymousName.str()] = mesh;
    return mesh;
  }
  
  // nElements was not specified, create and return anonymous standard regular mesh with 1 node, don't store it 
  std::array<element_idx_t, 1> nElements {0};
  std::array<double, 1> physicalExtent {1.0};
  LOG(DEBUG) << "Create new 1-node mesh with type "<<typeid(Mesh::RegularFixed<1>).name()<<", not stored.";
  std::shared_ptr<Mesh::Mesh> mesh = std::make_shared<Mesh::RegularFixed<1>>(nElements, physicalExtent);
  return mesh;
}
