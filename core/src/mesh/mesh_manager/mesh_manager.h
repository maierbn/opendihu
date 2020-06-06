#pragma once

#include <Python.h>  // has to be the first included header
#include <map>

#include "control/dihu_context.h"
//#include "function_space/function_space.h"
#include "mesh/mapping_between_meshes/manager/02_manager.h"
#include "function_space/function_space_generic.h"

namespace Partition{
class Manager;
}
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
  Manager(PythonConfig specificSettings);

  //! store the pointer to the partition manager
  void setPartitionManager(std::shared_ptr<Partition::Manager> partitionManager);

  //! return previously created mesh or create on the fly, already call functionSpace->initialize()
  template<typename FunctionSpaceType=FunctionSpace::Generic>
  std::shared_ptr<FunctionSpaceType> functionSpace(PythonConfig settings);
  
  //! return previously created mesh
  template<typename FunctionSpaceType=FunctionSpace::Generic>
  std::shared_ptr<FunctionSpaceType> functionSpace(std::string meshName);
  
  //! check if a function space with the given name and type is stored
  template<typename FunctionSpaceType>
  bool hasFunctionSpaceOfType(std::string meshName, bool outputWarning=false);

  //! check if a function space with the specified name is stored, the type is not checked
  bool hasFunctionSpace(std::string meshName);

  //! create a mesh not from python config but directly by calling an appropriate construtor.
  //! With this e.g. meshes from node positions can be created.
  //! Initialize the functionSpace from settings.
  template<typename FunctionSpaceType, typename ...Args>
  std::shared_ptr<FunctionSpaceType> createFunctionSpace(std::string name, Args && ...args);

  //! create a generic function space without mesh representation with dimension nEntries, in serial
  std::shared_ptr<FunctionSpace::Generic> createGenericFunctionSpace(int nEntries, std::string name);

  //! create a generic function space without mesh representation with dimension nEntries
  //! @param meshPartition This is only used to get the number of ranks.
  template<typename FunctionSpaceType>
  std::shared_ptr<FunctionSpace::Generic> createGenericFunctionSpace(int nEntries, std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, std::string name);

  //! create a mesh not from python config but directly by calling an appropriate construtor.
  //! With this e.g. meshes from node positions can be created.
  //! Use the given meshPartition, i.e. not the number of elements from settings.
  template<typename FunctionSpaceType, typename ...Args>
  std::shared_ptr<FunctionSpaceType> createFunctionSpaceWithGivenMeshPartition(
    std::string name, std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition, Args && ...args);

  //! Create a field variable without logical mesh representation, e.g. for MOR reduced vectors.
  //! The vector contains nEntries entries, the partitioning is done by the partition manager.
  //! \param nRanks on how many ranks the field variable will be distributed, the nEntries is the local size
  //! \param name is the name of the Petsc Vec, used for debugging output.
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::Generic,1>> createGenericFieldVariable(int nEntries, std::string name);

  friend class NodePositionsTester;    //< a class used for testing

private:

  //! Helper method to create a composite function space from the given settings, this uses the ManagerCompositeMesh helper class, in turn
  template<typename FunctionSpaceType>
  std::shared_ptr<FunctionSpaceType> createCompositeMesh(PythonConfig settings);

  struct NodePositionsFromFile
  {
    std::string filename;                            //< filename of the file to read
    std::vector<std::pair<MPI_Offset,int>> chunks;   //< pairs of (offset, number of values), where each value corresponds to 3 double values (position x,y,z) in data
    std::vector<double> data;                        //< the values of the node positions
  };

  //! store settings for all meshes that are specified in specificSettings_
  void storePreconfiguredMeshes();

  //! resolves the requested geometry data in nodePositionsFromFile_
  void loadGeometryFromFile();

  std::shared_ptr<Partition::Manager> partitionManager_;                //< the partition manager object
  PythonConfig specificSettings_;                                       //< the top level python settings
  
  int numberAnonymousMeshes_;                                           //< how many meshes without a given name in the python config are contained in meshes_. These have a key "anonymous<no>"

  std::map<std::string, PythonConfig> meshConfiguration_;               //< the python dicts for the meshes that were defined under "Meshes"
  std::map<std::string, std::shared_ptr<Mesh>> functionSpaces_;         //< the managed function spaces with their string key
  std::map<std::string, NodePositionsFromFile> nodePositionsFromFile_;  //< filename, offset, length, data of nodePosition data specified in a binary file
};

/** Helper class to create the composite meshes
 */
template<typename FunctionSpaceType>
class ManagerCompositeMesh
{
public:

  //! helper function to create a composite mesh
  static std::shared_ptr<FunctionSpaceType> createCompositeMesh(std::shared_ptr<Partition::Manager> partitionManager,
                                                                std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<FunctionSpaceType::dim()>,typename FunctionSpaceType::BasisFunction>>> subFunctionSpaces);
};

/** Helper class to create the composite meshes
 */
template<int D, typename BasisFunctionType>
class ManagerCompositeMesh<FunctionSpace::FunctionSpace<::Mesh::CompositeOfDimension<D>,BasisFunctionType>>
{
public:
  typedef FunctionSpace::FunctionSpace<::Mesh::CompositeOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! helper function to create a composite mesh
  static std::shared_ptr<FunctionSpaceType> createCompositeMesh(std::shared_ptr<Partition::Manager> partitionManager,
                                                                std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> subFunctionSpaces);
};

}  // namespace

#include "mesh/mesh_manager/mesh_manager.tpp"

