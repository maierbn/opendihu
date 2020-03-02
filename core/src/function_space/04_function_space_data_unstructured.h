#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/03_function_space_partition.h"
#include "mesh/type_traits.h"
#include "field_variable/unstructured/element_to_node_mapping.h"
#include "field_variable/00_field_variable_base.h"

// forward declaration of FieldVariable
namespace FieldVariable
{
template<typename FunctionSpaceType,int nComponents>
class FieldVariable;

template<typename FunctionSpaceType>
class FieldVariableBaseFunctionSpace;

template<typename FunctionSpaceType,int nComponents>
class Component;
}

namespace FunctionSpace
{
// forward declaration of FunctionSpace, needed for FieldVariable
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

/** Class that stores adjacency information for unstructured meshes and does file i/o for construction
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceDataUnstructured :
  public FunctionSpacePartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<FunctionSpaceDataUnstructured<D,BasisFunctionType>>
{
public:

  typedef ::FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> FieldVariableBaseFunctionSpaceType;
  typedef FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! constructor, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  FunctionSpaceDataUnstructured(std::shared_ptr<Partition::Manager> partitionManager, PythonConfig settings, bool noGeometryField=false);

  //! constructor, null argument is ignored
  FunctionSpaceDataUnstructured(std::shared_ptr<Partition::Manager> partitionManager, std::vector<double> &null, PythonConfig settings, bool noGeometryField=false);

  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! return the global/natural node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements, this is the same as getNodeNo because there are no global numbers for unstructured meshes
  global_no_t getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const;

  //! get all dofs of a specific node, as vector, the array version that is present for structured meshes does not make sense here, because with versions the number of dofs per node is not static.
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get the dof no of the first dof at the 
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;
  
  //! get the elementToNodeMapping
  std::shared_ptr<FieldVariable::ElementToNodeMapping> elementToNodeMapping();
  
  //! get the total number of elements on the local partition, for structured meshes this is directly implemented in the Mesh itself (not FunctionSpace like here)
  element_no_t nElementsLocal() const;

  //! get the total number of elements on the global domain, for structured meshes this is directly implemented in the Mesh itself (not FunctionSpace like here)
  global_no_t nElementsGlobal() const;

  //! get the node no in the global natural ordering
  global_no_t getNodeNoGlobalNaturalFromElementNoLocal(element_no_t elementNoLocal, int nodeIndex) const;

  // nDofsGlobal() is defined in 06_function_space_dofs_nodes.h
  
  //! initialize geometry
  virtual void initialize();
  
protected:
  //! parse a given *.exelem file and prepare fieldVariable_
  void parseExelemFile(std::string exelemFilename);

  //! parse a given *.exelem file and fill fieldVariable_ with values
  void parseExnodeFile(std::string exnodeFilename);

  //! rename field variables if "remap" is specified in config
  void remapFieldVariables(PythonConfig settings);

  //! multiply dof values with scale factors such that scale factor information is completely contained in dof values
  void eliminateScaleFactors();

  //! parse the element and node positions from python settings
  void parseFromSettings(PythonConfig settings);

  //! initialize the meshPartition of this mesh (by calling FunctionSpacePartition::initialize()), then create the partitioned Petsc vectors in each field variable
  void initializeValuesVector();
  
  std::map<std::string, std::shared_ptr<FieldVariableBaseFunctionSpaceType>> fieldVariable_; ///< all non-geometry field field variables that were present in exelem/exnode files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField_ = nullptr;  ///< the geometry field variable

  std::shared_ptr<FieldVariable::ElementToNodeMapping> elementToNodeMapping_;   ///< for every element the adjacent nodes and the field variable + dofs for their position
  element_no_t nElements_ = 0;    ///< number of elements in exelem file
  dof_no_t nDofs_ = 0;        ///< number of degrees of freedom. This can be different from nNodes * nDofsPerNode because of versions and shared nodes
  bool noGeometryField_;     ///< this is set if there is no geometry field stored. this is only needed for solid mechanics mixed formulation where the lower order basisOnMesh does not need its own geometry information

};

}  // namespace

#include "function_space/04_function_space_data_unstructured.tpp"
#include "function_space/04_function_space_data_unstructured_parse_exfiles.tpp"
#include "function_space/04_function_space_data_unstructured_parse_settings.tpp"
