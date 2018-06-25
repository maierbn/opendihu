#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "mesh/type_traits.h"
#include "mesh/face_t.h"

namespace BasisOnMesh
{

/** Base class to compute local dof and node no.s.
 * This is only possible for structured meshes because unstructured meshes store the information about
 * numbering explicitly.
 * The numberings refer to the stored information on this rank and is, thus, a local numbering.
 * It does not distinguish between ghost and interior nodes, also the ghost nodes are treated the same as every other stored node.
 */
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits=MeshType>
class BasisOnMeshNumbers
{};


/** class to compute local dof and node no.s, partial specialization for structured mesh, D=1
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> :
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;

  //! return the local dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
  
  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;
};

/** class to compute local dof and node no.s, partial specialization for structured mesh, D=2
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> :
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;

  //! return the local dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;

  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;
};

/** class to compute local dof and node no.s, partial specialization for structured mesh, D=3
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> :
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshJacobian<MeshType,BasisFunctionType>::BasisOnMeshJacobian;

  //! return the local dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;

  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;
};


}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_numbers_structured.tpp"