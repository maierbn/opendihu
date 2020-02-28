#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/04_function_space_numbers_structured.h"

namespace FunctionSpace
{

/** class to compute local dof and node no.s, partial specialization for composite mesh
 */

template<int D,typename BasisFunctionType>
class FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType> :
  public FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::FunctionSpacePartition;

  //! return the local dof number of element-local dof dofIndex of element elementNoLocal
  dof_no_t getDofNo(element_no_t elementNoLocal, int dofIndex) const;

  //! return the local node number of element-local node nodeIndex of element with local no elementNoLocal
  node_no_t getNodeNo(element_no_t elementNoLocal, int nodeIndex) const;

  //! get all dofs of a specific node, as vector
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;

  //! get all dofs of a specific node, as array
  void getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,FunctionSpaceBaseDim<D,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const;
  
  //! get the dof no of the specified dof at the node
  dof_no_t getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const;

  //! get neighbouring node to nodeNoLocal or -1 if there is no such node, nodeNoLocal has to be a non-ghost local node
  node_no_t getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction) const;

};

/** partial specialization for CompositeOfDimension<D>
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceNumbersCommon<Mesh::CompositeOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpaceNumbers<Mesh::CompositeOfDimension<D>,BasisFunctionType>::FunctionSpaceNumbers;

  //! get the node no in the global natural ordering
  global_no_t getNodeNoGlobalNaturalFromElementNoLocal(element_no_t elementNoLocal, int nodeIndex) const;

};

}  // namespace

#include "function_space/04_function_space_numbers_composite.tpp"
