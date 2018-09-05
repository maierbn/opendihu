#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/03_function_space_partition.h"
#include "mesh/type_traits.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

/** Base class to compute local dof and node no.s.
 * This is only possible for structured meshes because unstructured meshes store the information about
 * numbering explicitly.
 * The numberings refer to the stored information on this rank and is, thus, a local numbering.
 * It does not distinguish between ghost and interior nodes, also the ghost nodes are treated the same as every other stored node.
 */
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits=MeshType>
class FunctionSpaceNumbers
{};


/** class to compute local dof and node no.s, partial specialization for structured mesh, D=1
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> :
  public FunctionSpacePartition<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpacePartition<MeshType,BasisFunctionType>::FunctionSpacePartition;

  //! return the local dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element with local no elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! return the global/natural node number of element-local node nodeIndex of element with global no elementNoGlobal
  global_no_t getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
  
  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;

  //! get neighbouring node to nodeNoLocal or -1 if there is no such node, nodeNoLocal has to be a non-ghost local node
  node_no_t getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const;

  //! get node local no from the local coordinate in natural local numbering
  node_no_t getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const;
};

/** class to compute local dof and node no.s, partial specialization for structured mesh, D=2
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> :
  public FunctionSpacePartition<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpacePartition<MeshType,BasisFunctionType>::FunctionSpacePartition;

  //! return the local dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! return the global/natural node number of element-local node nodeIndex of element with global no elementNoGlobal
  global_no_t getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;

  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;

  //! get neighbouring node to nodeNoLocal or -1 if there is no such node, nodeNoLocal has to be a non-ghost local node
  node_no_t getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const;

  //! get node local no from the local coordinate in natural local numbering
  node_no_t getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const;
};

/** class to compute local dof and node no.s, partial specialization for structured mesh, D=3
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> :
  public FunctionSpacePartition<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpacePartition<MeshType,BasisFunctionType>::FunctionSpacePartition;

  //! return the local dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;

  //! return the global/natural node number of element-local node nodeIndex of element with global no elementNoGlobal
  global_no_t getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;

  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;

  //! get neighbouring node to nodeNoLocal or -1 if there is no such node, nodeNoLocal has to be a non-ghost local node
  node_no_t getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const;

  //! get node local no from the local coordinate in natural local numbering
  node_no_t getNodeNo(std::array<int,MeshType::dim()> coordinateLocal) const;
};


}  // namespace

#include "function_space/04_function_space_numbers_structured.tpp"
