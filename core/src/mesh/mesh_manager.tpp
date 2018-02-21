#include "mesh/mesh_manager.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/structured_regular_fixed.h"

namespace Mesh 
{

//! return previously created mesh or create on the fly
template<typename BasisOnMeshType>
std::shared_ptr<Mesh> Manager::mesh(PyObject *settings)
{     
  // if mesh was already created earlier
  if (PythonUtility::containsKey(settings, "meshName"))
  {
    std::string meshName = PythonUtility::getOptionString(settings, "meshName", "");
    if (hasMesh(meshName))
    {
      LOG(DEBUG) << "Mesh with meshName \""<<meshName<<"\" requested and found, type is "<<typeid(meshes_[meshName]).name();
      return std::static_pointer_cast<BasisOnMeshType>(meshes_[meshName]);
    }
    else if(meshConfiguration_.find(meshName) != meshConfiguration_.end())
    {
      // mesh was preconfigured, do nothing specific here, created standard mesh with 1 node
      LOG(DEBUG) << "Mesh configuration for \""<<meshName<<"\" found and requested, will be created now. "
        << "Type is "<< typeid(BasisOnMeshType).name()<<".";
      std::shared_ptr<BasisOnMeshType> mesh = std::make_shared<BasisOnMeshType>(meshConfiguration_[meshName]);
      mesh->initialize();
      meshes_[meshName] = mesh;
      LOG(DEBUG) << "Stored under key \""<<meshName<<"\".";
      return std::static_pointer_cast<BasisOnMeshType>(meshes_[meshName]);
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
  LOG(DEBUG) << "Create new mesh with type "<<typeid(BasisOnMeshType).name()<<" and name \""<<anonymousName.str()<<"\".";
  std::shared_ptr<BasisOnMeshType> mesh = std::make_shared<BasisOnMeshType>(settings);

  mesh->initialize();

  meshes_[anonymousName.str()] = mesh;
  
  VLOG(1) << "mesh nNodes: " << mesh->nNodes();
  
  return mesh;
}

}   // namespace