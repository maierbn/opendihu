#include "mesh/mesh_manager/mesh_manager.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/structured_regular_fixed.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace Mesh
{

//! return previously created mesh or create on the fly
template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> Manager::functionSpace(PythonConfig settings)
{
  LOG(DEBUG) << "querying Mesh::Manager::functionSpace, type " << StringUtility::demangle(typeid(FunctionSpaceType).name());
  
  // if mesh was already created earlier
  if (settings.hasKey("meshName"))
  {
    PyObject *meshNamePy = settings.getOptionPyObject("meshName");

    // if the value to "meshName" has multiple mesh names, this is to indicate a composite mesh
    if (PyList_Check(meshNamePy))
    {
      return this->createCompositeMesh<FunctionSpaceType>(settings);
    }
    else
    {
      // this is not a composite mesh, but just a normal mesh, get the meshName and create the mesh
      std::string meshName = settings.getOptionString("meshName", "");

      // output error message if for a composite mesh the meshName is no list
      if (isComposite<std::shared_ptr<FunctionSpaceType>>::value)
      {
        LOG(FATAL) << "Function space of type " << StringUtility::demangle(typeid(FunctionSpaceType).name())
          << " is composite and requires a list of mesh names under " << settings << "[\"meshName\"]. "
          << "Currently, \"meshName\" is not a list, but the entry \"" << meshName << "\".";
      }
      std::shared_ptr<FunctionSpaceType> functionSpace = this->functionSpace<FunctionSpaceType>(meshName);

      if (functionSpace)
      {
        return functionSpace;
      }
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
  LOG(DEBUG) << "Create new mesh with type " << StringUtility::demangle(typeid(FunctionSpaceType).name()) << " and name \"" <<anonymousName.str() << "\" (1).";
  
  // create mesh and initialize
  std::shared_ptr<FunctionSpaceType> functionSpace = createFunctionSpace<FunctionSpaceType>(anonymousName.str(), settings);
  
  std::string logKey;
  if (settings.hasKey("logKey"))
  {
    logKey = settings.getOptionString("logKey", "");
  }
  
  Control::PerformanceMeasurement::setParameter(std::string("nDofs") + logKey, functionSpace->nDofsGlobal());
  Control::PerformanceMeasurement::setParameter(std::string("nNodes") + logKey, functionSpace->nNodesGlobal());
  Control::PerformanceMeasurement::setParameter(std::string("nElements") + logKey, functionSpace->nElementsGlobal());
  
  return functionSpace;
}

//! return previously created mesh or create on the fly
template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> Manager::functionSpace(std::string meshName)
{
  LOG(DEBUG) << "querying Mesh::Manager::functionSpace, type " << StringUtility::demangle(typeid(FunctionSpaceType).name());

  if (hasFunctionSpaceOfType<FunctionSpaceType>(meshName))
  {
    LOG(DEBUG) << "Mesh with meshName \"" << meshName << "\" requested, found and type matches, type is "
      << StringUtility::demangle(typeid(functionSpaces_[meshName]).name())
    << ", cast to " << StringUtility::demangle(typeid(FunctionSpaceType).name());
    return std::static_pointer_cast<FunctionSpaceType>(functionSpaces_[meshName]);
  }
  else if (meshConfiguration_.find(meshName) != meshConfiguration_.end())
  {
    // mesh was preconfigured, create new mesh from stored meshConfiguration
    LOG(DEBUG) << "Mesh configuration for \"" << meshName << "\" found and requested, will be created now. "
    << "Type is " << StringUtility::demangle(typeid(FunctionSpaceType).name()) << ".";

    // get mesh configuration that was parsed earlier
    PythonConfig meshConfiguration = meshConfiguration_.at(meshName);

    if (hasFunctionSpace(meshName))
    {
      std::stringstream newMeshName;
      newMeshName << meshName << "_2";
      meshName = newMeshName.str();
      LOG(INFO) << "Create a mesh with name \"" << meshName << "\".";
    }

    // create new mesh and initialize
    std::shared_ptr<FunctionSpaceType> functionSpace = createFunctionSpace<FunctionSpaceType>(meshName, meshConfiguration);

    std::string logKey;
    if (meshConfiguration.hasKey("logKey"))
    {
      logKey = meshConfiguration.getOptionString("logKey", "");
    }

    Control::PerformanceMeasurement::setParameter(std::string("~nDofs") + logKey, functionSpace->nDofsGlobal());
    Control::PerformanceMeasurement::setParameter(std::string("~nNodes") + logKey, functionSpace->nNodesGlobal());
    Control::PerformanceMeasurement::setParameter(std::string("~nElements") + logKey, functionSpace->nElementsGlobal());

    return functionSpace;
  }
  else
  {
    LOG(ERROR) << "Mesh with name \"" << meshName << "\" was not defined. "
      << "Add the following to the python settings:\n\"Meshes\": {\n\t\"" << meshName << "\": { \n\t\t<your mesh options here>\n\t}\n}";
  }

  return nullptr;
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
  LOG(DEBUG) << "Create new mesh with type " << StringUtility::demangle(typeid(FunctionSpaceType).name())
    << " and name \"" <<name << "\" (2).";

  // create mesh and initialize
  std::shared_ptr<FunctionSpaceType> functionSpace;

  // check if node positions from file are available
  if (nodePositionsFromFile_.find(name) != nodePositionsFromFile_.end())
  {
    std::vector<double> &nodePositions = nodePositionsFromFile_[name].data;
    functionSpace = std::make_shared<FunctionSpaceType>(this->partitionManager_, nodePositions, std::forward<Args>(args)...);
  }
  else
  {
    // construct normally
    functionSpace = std::make_shared<FunctionSpaceType>(this->partitionManager_, std::forward<Args>(args)...);
  }

  functionSpace->setMeshName(name);
  functionSpace->initialize();

  functionSpaces_[name] = functionSpace;

  VLOG(1) << "mesh nNodes: (local without ghosts: " << functionSpace->nNodesLocalWithoutGhosts()
    << ", with ghosts: " << functionSpace->nNodesLocalWithGhosts() << "), stored under key \"" << name << "\"";

  return functionSpace;
}

template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> Manager::createCompositeMesh(PythonConfig settings)
{
  std::vector<std::string> meshNames;
  settings.template getOptionVector<std::string>("meshName", meshNames);

  LOG(DEBUG) << "create composite Mesh from meshes " << meshNames;

  // create the meshName of the composite mesh as concatenation of all sub mesh names with the '+' sign
  std::stringstream compositeName;
  int i = 0;
  for (std::vector<std::string>::iterator iter = meshNames.begin(); iter != meshNames.end(); iter++, i++)
  {
    if (i != 0)
      compositeName << "+";
    compositeName << *iter;
  }

  // check if this mesh has already been created
  if (hasFunctionSpaceOfType<FunctionSpaceType>(compositeName.str()))
  {
    LOG(DEBUG) << "Composite mesh with meshName \"" << compositeName.str() << "\" requested, found and type matches, type is "
      << StringUtility::demangle(typeid(functionSpaces_[compositeName.str()]).name())
    << ", cast to " << StringUtility::demangle(typeid(FunctionSpaceType).name());
    return std::static_pointer_cast<FunctionSpaceType>(functionSpaces_[compositeName.str()]);
  }
  else
  {
    // create sub meshes
    typedef FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<FunctionSpaceType::dim()>,typename FunctionSpaceType::BasisFunction> SubFunctionSpaceType;
    std::vector<std::shared_ptr<SubFunctionSpaceType>> subFunctionSpaces;
    int i = 0;
    for (std::vector<std::string>::iterator iter = meshNames.begin(); iter != meshNames.end(); iter++, i++)
    {
      LOG(DEBUG) << "For composite mesh create submesh " << i << "/" << meshNames.size();
      subFunctionSpaces.push_back(this->functionSpace<SubFunctionSpaceType>(*iter));
    }

    // create new composite mesh and initialize
    std::shared_ptr<FunctionSpaceType> functionSpace = ManagerCompositeMesh<FunctionSpaceType>::createCompositeMesh(this->partitionManager_, subFunctionSpaces);

    functionSpace->setMeshName(compositeName.str());
    functionSpace->initialize();

    // store created mesh
    functionSpaces_[compositeName.str()] = functionSpace;

    // add statistics information to log
    std::string logKey = compositeName.str();
    Control::PerformanceMeasurement::setParameter(std::string("~nDofs") + logKey, functionSpace->nDofsGlobal());
    Control::PerformanceMeasurement::setParameter(std::string("~nNodes") + logKey, functionSpace->nNodesGlobal());
    Control::PerformanceMeasurement::setParameter(std::string("~nElements") + logKey, functionSpace->nElementsGlobal());

    return functionSpace;
  }
}

template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> ManagerCompositeMesh<FunctionSpaceType>::
createCompositeMesh(std::shared_ptr<Partition::Manager> partitionManager,
                    std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<FunctionSpaceType::dim()>,typename FunctionSpaceType::BasisFunction>>> subFunctionSpaces)
{
  LOG(FATAL) << "Trying to create a composite function space from the Python Confing but the C++ mesh type is not Mesh::CompositeOfDimension<D> "
    << "but " << StringUtility::demangle(typeid(typename FunctionSpaceType::Mesh).name());
  return nullptr;
}

template<int D, typename BasisFunctionType>
std::shared_ptr<FunctionSpace::FunctionSpace<::Mesh::CompositeOfDimension<D>,BasisFunctionType>>
ManagerCompositeMesh<FunctionSpace::FunctionSpace<::Mesh::CompositeOfDimension<D>,BasisFunctionType>>::
createCompositeMesh(std::shared_ptr<Partition::Manager> partitionManager,
                    std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> subFunctionSpaces)
{
  // create the composite function space, directly using the subFunctionSpaces
  return std::make_shared<FunctionSpace::FunctionSpace<::Mesh::CompositeOfDimension<D>,BasisFunctionType>>(partitionManager, subFunctionSpaces);
}

template<typename FunctionSpaceType, typename ...Args>
std::shared_ptr<FunctionSpaceType> Manager::createFunctionSpaceWithGivenMeshPartition(
  std::string name, std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, Args && ...args)
{
  if (hasFunctionSpaceOfType<FunctionSpaceType>(name))
  {
    LOG(ERROR) << "FunctionSpace with name \"" <<name << "\" already exists. Overwrite mesh.";
  }

  // create new mesh
  LOG(DEBUG) << "Create new mesh with type " << StringUtility::demangle(typeid(FunctionSpaceType).name())
    << " and name \"" <<name << "\" (3).";

  // create mesh and initialize
  std::shared_ptr<FunctionSpaceType> functionSpace = std::make_shared<FunctionSpaceType>(this->partitionManager_, std::forward<Args>(args)...);

  functionSpace->setMeshName(name);
  functionSpace->setMeshPartition(meshPartition);
  functionSpace->initialize();

  functionSpaces_[name] = functionSpace;

  VLOG(1) << "mesh nNodes: (local without ghosts: " << functionSpace->nNodesLocalWithoutGhosts()
    << ", with ghosts: " << functionSpace->nNodesLocalWithGhosts() << "), stored under key \"" << name << "\"";

  return functionSpace;
}

template<typename FunctionSpaceType>
bool Manager::hasFunctionSpaceOfType(std::string meshName)
{
  VLOG(1) << "hasMesh(" << meshName << "): " << (functionSpaces_.find(meshName) != functionSpaces_.end());
  VLOG(1) << "meshes size: " << functionSpaces_.size();

  if (functionSpaces_.find(meshName) != functionSpaces_.end()) // if mesh is found by name
  {
    if (std::dynamic_pointer_cast<FunctionSpaceType>(functionSpaces_[meshName]))
    {
      return true;
    }
    else
    {
      LOG(DEBUG) << "Mesh \"" << meshName << "\" is stored but under a type that is different from "
        << StringUtility::demangle(typeid(FunctionSpaceType).name()) << ". A new mesh will be created instead.";
      // This warning could be if in an operator splitting setup for Cellml adapter the functionType is set differently from the one in the FEM.
    }
  }

  return false;
}

template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpace::Generic> Manager::
createGenericFunctionSpace(int nEntries, std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, std::string name)
{
  std::array<element_no_t, 1> nElements;

  if (meshPartition->ownRankNo() == meshPartition->nRanks()-1)
  {
    // on the right border, have one element less than nodes, because nodes are left and right of the 1D element
    nElements[0] = nEntries-1;
  }
  else
  {
    // everywhere else there are as many nodes as elements
    nElements[0] = nEntries;
  }

  LOG(DEBUG) << "createGenericFunctionSpace, nEntries: " << nEntries;
  // << ", hasFullNumberOfNodes: "
  //  << "["  << meshPartition->hasFullNumberOfNodes(0) << ", "  << meshPartition->hasFullNumberOfNodes(1) << ", "  << meshPartition->hasFullNumberOfNodes(2) << "], nElements: " << nElements;

  std::array<double, 1> physicalExtent({0.0});
  std::array<int, 1> nRanksPerCoordinateDirection({meshPartition->nRanks()});
  std::shared_ptr<Mesh> mesh = createFunctionSpace<FunctionSpace::Generic>(name, nElements, physicalExtent, nRanksPerCoordinateDirection, false);   // last parameter is that nElements is local number

  return std::static_pointer_cast<FunctionSpace::Generic>(mesh);
}
}   // namespace
