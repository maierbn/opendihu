#include "mesh/mesh_manager.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/structured_regular_fixed.h"

namespace Mesh
{

//! return previously created mesh or create on the fly
template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> Manager::functionSpace(PyObject *settings)
{
  LOG(DEBUG) << "querying Mesh::Manager::functionSpace, type " << typeid(FunctionSpaceType).name();

  // if mesh was already created earlier
  if (PythonUtility::hasKey(settings, "meshName"))
  {
    std::string meshName = PythonUtility::getOptionString(settings, "meshName", "");
    if (hasFunctionSpaceOfType<FunctionSpaceType>(meshName))
    {
      LOG(DEBUG) << "Mesh with meshName \"" << meshName << "\" requested, found and type matches, type is " << typeid(functionSpaces_[meshName]).name()
        << ", cast to " << typeid(FunctionSpaceType).name();
      return std::static_pointer_cast<FunctionSpaceType>(functionSpaces_[meshName]);
    }
    else if(meshConfiguration_.find(meshName) != meshConfiguration_.end())
    {
      // mesh was preconfigured, create new mesh from stored meshConfiguration
      LOG(DEBUG) << "Mesh configuration for \"" << meshName << "\" found and requested, will be created now. "
        << "Type is " << typeid(FunctionSpaceType).name() << ".";
      
      // get mesh configuration that was parsed earlier
      PyObject *meshConfiguration = meshConfiguration_.at(meshName);
      
      if (hasFunctionSpace(meshName))
      {
        std::stringstream newMeshName;
        newMeshName << meshName << "_2";
        meshName = newMeshName.str();
        LOG(INFO) << "Create a mesh with name \"" << meshName << "\".";
      }

      // create new mesh and initialize
      std::shared_ptr<FunctionSpaceType> functionSpace = std::make_shared<FunctionSpaceType>(this->partitionManager_, meshConfiguration);
      functionSpace->setMeshName(meshName);
      functionSpace->initialize();
      
      // store mesh under its name
      functionSpaces_[meshName] = functionSpace;
      LOG(DEBUG) << "Stored under key \"" << meshName << "\".";
      return functionSpace;
    }
    else
    {
      LOG(ERROR) << "Config contains reference to mesh with meshName \"" << meshName << "\" but no such mesh was defined.";
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
  LOG(DEBUG) << "Create new mesh with type " << typeid(FunctionSpaceType).name() << " and name \"" <<anonymousName.str() << "\".";
  
  // create mesh and initialize
  std::shared_ptr<FunctionSpaceType> functionSpace = std::make_shared<FunctionSpaceType>(this->partitionManager_, settings);
  functionSpace->setMeshName(anonymousName.str());
  functionSpace->initialize();

  functionSpaces_[anonymousName.str()] = functionSpace;

  VLOG(1) << "mesh nNodes: (local without ghosts: " << functionSpace->nNodesLocalWithoutGhosts() << ", with ghosts: " << functionSpace->nNodesLocalWithGhosts() << ")";

  return functionSpace;
}

//! create a mesh not from python config but directly by calling an appropriate construtor. 
//! With this e.g. meshes from node positions can be created.
template<typename FunctionSpaceType, typename ...Args>
std::shared_ptr<FunctionSpaceType> Manager::createFunctionSpace(std::string name, Args && ...args)
{
  if (hasFunctionSpaceOfType<FunctionSpaceType>(name))
  {
    LOG(ERROR) << "FunctionSpace with name \"" <<name << "\" already exists. Overwrite mesh.";
  }
 
  // create new mesh
  LOG(DEBUG) << "Create new mesh with type " << typeid(FunctionSpaceType).name() << " and name \"" <<name << "\".";
  
  // create mesh and initialize
  std::shared_ptr<FunctionSpaceType> functionSpace = std::make_shared<FunctionSpaceType>(this->partitionManager_, std::forward<Args>(args)...);

  functionSpace->setMeshName(name);
  functionSpace->initialize();
  
  functionSpaces_[name] = functionSpace;

  VLOG(1) << "mesh nNodes: (local without ghosts: " << functionSpace->nNodesLocalWithoutGhosts()
    << ", with ghosts: " << functionSpace->nNodesLocalWithGhosts() << ")";

  return functionSpace;
}

template<typename FunctionSpaceType>
bool Manager::hasFunctionSpaceOfType(std::string meshName)
{
  LOG(DEBUG) << "hasMesh(" << meshName << "): " << (functionSpaces_.find(meshName) != functionSpaces_.end());
  LOG(DEBUG) << "meshes size: " << functionSpaces_.size();

  if (functionSpaces_.find(meshName) != functionSpaces_.end()) // if mesh is found by name
  {
    if (std::dynamic_pointer_cast<FunctionSpaceType>(functionSpaces_[meshName]))
    {
      return true;
    }
    else
    {
      LOG(WARNING) << "Mesh \"" << meshName << "\" is stored but under a type that is different from " << typeid(FunctionSpaceType).name() << ". A new mesh will be created instead.";
      // This warning could be if in an operator splitting setup for Cellml adapter the functionType is set differently from the one in the FEM.
    }
  }

  return false;
}

}   // namespace
