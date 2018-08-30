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
  if (PythonUtility::hasKey(settings, "meshName"))
  {
    std::string meshName = PythonUtility::getOptionString(settings, "meshName", "");
    if (hasMesh(meshName))
    {
      LOG(DEBUG) << "Mesh with meshName \"" <<meshName << "\" requested and found, type is " << typeid(meshes_[meshName]).name();
      return std::static_pointer_cast<BasisOnMeshType>(meshes_[meshName]);
    }
    else if(meshConfiguration_.find(meshName) != meshConfiguration_.end())
    {
      // mesh was preconfigured, create new mesh from stored meshConfiguration
      LOG(DEBUG) << "Mesh configuration for \"" <<meshName << "\" found and requested, will be created now. "
        << "Type is " << typeid(BasisOnMeshType).name() << ".";
      
      // get mesh configuration that was parsed earlier
      PyObject *meshConfiguration = meshConfiguration_.at(meshName);
      
      // create new mesh and initialize
      std::shared_ptr<BasisOnMeshType> mesh = std::make_shared<BasisOnMeshType>(this->partitionManager_, meshConfiguration);
      mesh->setMeshName(meshName);
      mesh->initialize();
      
      // store mesh under its name
      meshes_[meshName] = mesh;
      LOG(DEBUG) << "Stored under key \"" <<meshName << "\".";
      return std::static_pointer_cast<BasisOnMeshType>(meshes_[meshName]);
    }
    else
    {
      LOG(ERROR) << "Config contains reference to mesh with meshName \"" <<meshName << "\" but no such mesh was defined.";
    }
  }
  else
  {
    LOG(DEBUG) << "Config does not contain meshName.";
  }

  // if there was no name given for the mesh
  
  // create new mesh, store as anonymous object
  std::stringstream anonymousName;
  anonymousName << "anonymous" << numberAnonymousMeshes_++;
  LOG(DEBUG) << "Create new mesh with type " << typeid(BasisOnMeshType).name() << " and name \"" <<anonymousName.str() << "\".";
  
  // create mesh and initialize
  std::shared_ptr<BasisOnMeshType> mesh = std::make_shared<BasisOnMeshType>(this->partitionManager_, settings);
  mesh->setMeshName(anonymousName.str());
  mesh->initialize();

  meshes_[anonymousName.str()] = mesh;

  VLOG(1) << "mesh nNodes: (local without ghosts: " << mesh->nNodesLocalWithoutGhosts() << ", with ghosts: " << mesh->nNodesLocalWithGhosts() << ")";

  return mesh;
}

//! create a mesh not from python config but directly by calling an appropriate construtor. 
//! With this e.g. meshes from node positions can be created.
template<typename BasisOnMeshType, typename ...Args>
std::shared_ptr<Mesh> Manager::createMesh(std::string name, Args && ...args)
{
  if (this->hasMesh(name))
  {
    LOG(ERROR) << "Mesh with name \"" <<name << "\" already exists. Overwrite mesh.";
  }
 
  // create new mesh
  LOG(DEBUG) << "Create new mesh with type " << typeid(BasisOnMeshType).name() << " and name \"" <<name << "\".";
  
  // create mesh and initialize
  std::shared_ptr<BasisOnMeshType> mesh = std::make_shared<BasisOnMeshType>(this->partitionManager_, std::forward<Args>(args)...);

  mesh->setMeshName(name);
  mesh->initialize();
  
  meshes_[name] = mesh;

  VLOG(1) << "mesh nNodes: (local without ghosts: " << mesh->nNodesLocalWithoutGhosts() << ", with ghosts: " << mesh->nNodesLocalWithGhosts() << ")";

  return mesh;
}

}   // namespace