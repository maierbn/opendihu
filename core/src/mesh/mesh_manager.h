#pragma once

#include <Python.h>  // has to be the first included header
#include <map>

#include "control/dihu_context.h"
#include "basis_on_mesh/basis_on_mesh.h"

namespace Mesh
{
class NodePositionsTester;

/**
 * This class creates and stores all used meshes.
 * Each mesh can be defined in the python config under "Meshes" with a name and other properties.
 * Various components of the program can later
 * request their mesh by a call to mesh(name).
 * If a mesh was not defined earlier, it is created on the fly when it is requested.
 */
class Manager
{
public:
  //! constructor
  Manager(DihuContext context);
  
  //! return previously created mesh or create on the fly
  template<typename BasisOnMeshType=None>
  std::shared_ptr<Mesh> mesh(PyObject *settings);
  
  //! check if a mesh with the given name is stored
  bool hasMesh(std::string meshName);
  
  friend class NodePositionsTester;    ///< a class used for testing 
  
private:
  //! store settings for all meshes that are specified in specificSettings_
  void storePreconfiguredMeshes();
 
  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key, for meshManager
  int numberAnonymousMeshes_;     ///< how many meshes without a given name in the python config are contained in meshes_. These have a key "anonymous<no>"
  
  std::map<std::string, PyObject *> meshConfiguration_;         ///< the python dicts for the meshes that were defined under "Meshes"
  std::map<std::string, std::shared_ptr<Mesh>> meshes_;    ///< the managed meshes with their string key
};

};    // namespace

#include "mesh/mesh_manager.tpp"
