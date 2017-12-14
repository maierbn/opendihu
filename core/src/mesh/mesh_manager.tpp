#include "mesh/mesh_manager.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/regular_fixed.h"

template<>
std::shared_ptr<Mesh::Mesh> MeshManager::mesh<Mesh::Mesh>(PyObject *settings);

//! return previously created mesh or create on the fly
template<class MeshType>
std::shared_ptr<MeshType> MeshManager::mesh(PyObject *settings)
{
  // if mesh was already created earlier
  if (PythonUtility::containsKey(settings, "meshName"))
  {
    std::string meshName = PythonUtility::getOptionString(settings, "meshName", "");
    if (hasMesh(meshName))
    {
      LOG(DEBUG) << "Mesh with meshName \""<<meshName<<"\" requested and found, type is "<<typeid(meshes_[meshName]).name();
      return std::static_pointer_cast<MeshType>(meshes_[meshName]);
    }
    else if(meshConfiguration_.find(meshName) != meshConfiguration_.end())
    {
      // mesh was preconfigured, do nothing specific here, created standard mesh with 1 node
      LOG(DEBUG) << "Mesh configuration for \""<<meshName<<"\" found and requested, will be created now. "
        << "Type is "<< typeid(MeshType).name()<<".";
      meshes_[meshName] = std::make_shared<MeshType>(meshConfiguration_[meshName]);
      LOG(DEBUG) << "Stored under key \""<<meshName<<"\".";
      return std::static_pointer_cast<MeshType>(meshes_[meshName]);
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
  
  // create new mesh, store as anonymous object
  std::stringstream anonymousName;
  anonymousName << "anonymous" << numberAnonymousMeshes_++;
  LOG(DEBUG) << "Create new mesh with type "<<typeid(MeshType).name()<<" and name \""<<anonymousName.str()<<"\".";
  std::shared_ptr<MeshType> mesh = std::make_shared<MeshType>(settings);
  meshes_[anonymousName.str()] = mesh;
  return mesh;
}