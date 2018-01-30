#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "mesh/type_traits.h"
#include "field_variable/element_to_node_mapping.h"
//#include "field_variable/field_variable.h"

// forward declaration of FieldVariable
namespace FieldVariable
{
template<typename BasisOnMeshType>
class FieldVariable;
};

namespace BasisOnMesh
{
// forward declaration of BasisOnMesh, needed for FieldVariable
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh;

/** base class to compute global dof and node no.s, for unstructured meshes
 */
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits=MeshType>
class BasisOnMeshDofs : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
  
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  int getDofNo(element_idx_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  int getNodeNo(element_idx_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const;
};

/** partial specialization for structured mesh, D=1
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
  
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the number of elements in each coordinate direction
  static int getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex);
  
  //! return the global dof number of element-local dof dofIndex of element elementNo
  int getDofNo(element_idx_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  int getNodeNo(element_idx_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const;
};

/** partial specialization for structured mesh, D=2
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
 
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the number of elements in each coordinate direction
  static int getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex);
  
  //! return the global dof number of element-local dof dofIndex of element elementNo
  int getDofNo(element_idx_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  int getNodeNo(element_idx_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const;
};

/** partial specialization for structured mesh, D=3
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
 
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the number of elements in each coordinate direction
  static int getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex);
  
  //! return the global dof number of element-local dof dofIndex of element elementNo
  int getDofNo(element_idx_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  int getNodeNo(element_idx_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const;
};

/** partial specialization for unstructured mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType> :
  public BasisOnMeshJacobian<Mesh::UnstructuredDeformable<D>,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>
{
public:

  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> FieldVariableType;
 
  //! constructor
  BasisOnMeshDofs(PyObject *settings);
  
  //! this assigns the geometry field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();
  
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  int getDofNo(element_idx_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  int getNodeNo(element_idx_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const;
  
  //! write exelem file to stream
  void outputExelemFile(std::ostream &file);
  
  //! write exnode file to stream
  void outputExnodeFile(std::ostream &file);
  
  //! return the number of elements
  element_idx_t nElements() const;
  
  //! return the number of dofs
  dof_idx_t nDofs() const;
 
  //! add all field variables except the geometry field to the vector. This is used to retrive additional field variables that were parsed from an exfile.
  void addNonGeometryFieldVariables(std::vector<std::shared_ptr<FieldVariableType>> &fieldVariables);
  
protected:
  //! parse a given *.exelem file and prepare fieldVariable_
  void parseExelemFile(std::string exelemFilename);
  
  //! parse a given *.exelem file and fill fieldVariable_ with values
  void parseExnodeFile(std::string exnodeFilename);
  
  //! rename field variables if "remap" is specified in config
  void remapFieldVariables(PyObject *settings);
  
  //! multiply dof values with scale factors such that scale factor information is completely contained in dof values
  void eliminateScaleFactors();
  
  //! parse the element and node positions from python settings
  void parseFromSettings(PyObject *settings);
  
  std::map<std::string, std::shared_ptr<FieldVariableType>> fieldVariable_;   ///< all field variables that were present in exelem/exnode files, should contain "geometry" field variable
  std::shared_ptr<FieldVariable::ElementToNodeMapping> elementToNodeMapping_;   ///< for every element the adjacent nodes and the field variable + dofs for their position
  element_idx_t nElements_ = 0;    ///< number of elements in exelem file 
  dof_idx_t nDofs_ = 0;        ///< number of degrees of freedom. This can be different from nNodes * nDofsPerNode because of versions and shared nodes
 
}; 
}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_dofs.tpp"
#include "basis_on_mesh/03_basis_on_mesh_dofs_input_output.tpp"
