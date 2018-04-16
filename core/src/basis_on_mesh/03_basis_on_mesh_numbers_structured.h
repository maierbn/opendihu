#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "mesh/type_traits.h"

namespace BasisOnMesh
{

/** base class to compute global dof and node no.s. 
 * This is only possible for structured meshes because unstructured meshes store the information about 
 * numbering explicitly.
 */
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits=MeshType>
class BasisOnMeshNumbers
{};


/** class to compute global dof and node no.s, partial specialization for structured mesh, D=1
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
  
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  
  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
};

/** class to compute global dof and node no.s, partial specialization for structured mesh, D=2
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
 
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  
  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
};

/** class to compute global dof and node no.s, partial specialization for structured mesh, D=3
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> : 
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;
 
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;  
  
  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
};


}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_numbers_structured.tpp"