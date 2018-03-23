#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "mesh/type_traits.h"
#include "field_variable/unstructured/element_to_node_mapping.h"
#include "field_variable/00_field_variable_base.h"

// forward declaration of FieldVariable
namespace FieldVariable
{
template<typename BasisOnMeshType,int nComponents>
class FieldVariable;

template<typename BasisOnMeshType>
class FieldVariableBase;

template<typename BasisOnMeshType>
class Component;
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
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  
  typedef ::FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>> FieldVariableBaseType;
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name){return nullptr;}
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
  
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
  
  typedef ::FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>> FieldVariableBaseType;
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name){return nullptr;}
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
 
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
  
  typedef ::FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>> FieldVariableBaseType;
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name){return nullptr;}
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
 
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;  
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
  
  typedef ::FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>> FieldVariableBaseType;
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name){return nullptr;}
};

/** partial specialization for unstructured mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshJacobian<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:

  typedef ::FieldVariable::FieldVariableBase<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> FieldVariableBaseType;
  typedef BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! constructor, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  BasisOnMeshDofs(PyObject *settings, bool noGeometryField=false);
  
  //! this assigns the geometry field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();
  
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! get all dofs of a single node  
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  
  //! write exelem file to stream
  void outputExelemFile(std::ostream &file);
  
  //! write exnode file to stream
  void outputExnodeFile(std::ostream &file);
  
  //! return the number of elements
  element_no_t nElements() const;
  
  //! return the number of dofs
  dof_no_t nDofs() const;
 
  //! access all stored field variables, the type FieldVariableBaseType does not contain the nComponents template paramet of FieldVariableType, so a cast should be performed in order to access the field variable
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name);
  
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
  
  std::map<std::string, std::shared_ptr<FieldVariableBaseType>> fieldVariable_; ///< all non-geometry field field variables that were present in exelem/exnode files
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> geometryField_ = nullptr;  ///< the geometry field variable
  
  std::shared_ptr<FieldVariable::ElementToNodeMapping> elementToNodeMapping_;   ///< for every element the adjacent nodes and the field variable + dofs for their position
  element_no_t nElements_ = 0;    ///< number of elements in exelem file 
  dof_no_t nDofs_ = 0;        ///< number of degrees of freedom. This can be different from nNodes * nDofsPerNode because of versions and shared nodes
  bool noGeometryField_;     ///< this is set if there is no geometry field stored. this is only needed for solid mechanics mixed formulation where the lower order basisOnMesh does not need its own geometry information
 
}; 

/** partial specialization for structured mesh and complete polynomials
 */
template<typename MeshType,int order>
class BasisOnMeshDofs<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>,Mesh::isDeformable<MeshType>> :
  public BasisOnMeshFunction<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:
  //! inherit constructor
  using BasisOnMeshFunction<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::BasisOnMeshJacobian;
  
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
};

}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_dofs_structured.tpp"
#include "basis_on_mesh/03_basis_on_mesh_dofs_unstructured.tpp"
#include "basis_on_mesh/03_basis_on_mesh_dofs_input_output.tpp"
