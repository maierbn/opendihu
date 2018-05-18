#include "mesh/mesh_manager.h"

#include "basis_on_mesh/basis_on_mesh.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/unstructured_deformable.h"

namespace Mesh
{

Manager::Manager(PyObject *specificSettings) :
  specificSettings_(specificSettings), numberAnonymousMeshes_(0)
{
  LOG(TRACE) << "MeshManager constructor";
  storePreconfiguredMeshes();
}

void Manager::storePreconfiguredMeshes()
{
  LOG(TRACE) << "MeshManager::storePreconfiguredMeshes";
  if (specificSettings_)
  {
    std::string keyString("Meshes");
    if (PythonUtility::hasKey(specificSettings_, "Meshes"))
    {

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
    else
    {
      LOG(INFO) << "You have specified the mesh in-line and not under the extra key \"Meshes\". You could do so,"
        " by defining \"Meshes\": {\"<your custom mesh name>\": {<your mesh parameters>}} at the beginning of the"
        " config and \"meshName\": \"<your custom mesh name>\" where you currently have specified the mesh parameters."
        " This is required if you want to use the same mesh for multiple objects.";
    }
  }
}

bool Manager::hasMesh(std::string meshName)
{
  LOG(DEBUG) << "hasMesh("<<meshName<<")";
  LOG(DEBUG) << "meshes size: " << meshes_.size();

  return meshes_.find(meshName) != meshes_.end();
}

template<>
std::shared_ptr<Mesh> Manager::
mesh<None>(PyObject *settings)
{
  std::string meshName;
  if (PythonUtility::hasKey(settings, "meshName"))
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
        << " Type is not clear, so go for StructuredRegularFixedOfDimension<1>.";
      typedef BasisOnMesh::BasisOnMesh<StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<>> NewBasisOnMesh;
      std::shared_ptr<NewBasisOnMesh> mesh = std::make_shared<NewBasisOnMesh>(meshConfiguration_[meshName]);
      mesh->initialize();
      meshes_[meshName] = mesh;
      LOG(DEBUG) << "Stored under key \""<<meshName<<"\".";
      return std::static_pointer_cast<Mesh>(meshes_[meshName]);
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
  if (PythonUtility::hasKey(settings, "nElements"))
  {
    // create new mesh, store as anonymous object
    std::stringstream anonymousName;
    anonymousName << "anonymous" << numberAnonymousMeshes_++;

    // set type to be 1D regular fixed mesh with linear lagrange basis
    typedef BasisOnMesh::BasisOnMesh<StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<>> NewBasisOnMesh;
    LOG(DEBUG) << "Create new mesh with type "<<typeid(NewBasisOnMesh).name()<<" and name \""<<anonymousName.str()<<"\".";

    std::shared_ptr<NewBasisOnMesh> mesh = std::make_shared<NewBasisOnMesh>(settings);
    mesh->initialize();

    meshes_[anonymousName.str()] = mesh;
    return mesh;
  }

  // nElements was not specified, create and return anonymous standard regular mesh with 1 node, don't store it
  std::array<element_no_t, 1> nElements {0};
  std::array<double, 1> physicalExtent {1.0};

  typedef BasisOnMesh::BasisOnMesh<StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<>> NewBasisOnMesh;
  LOG(DEBUG) << "Create new 1-node mesh with type "<<typeid(NewBasisOnMesh).name()<<", not stored.";

  std::shared_ptr<NewBasisOnMesh> mesh = std::make_shared<NewBasisOnMesh>(nElements, physicalExtent);
  mesh->initialize();
  return mesh;
}

};  // namespace