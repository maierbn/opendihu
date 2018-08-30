#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/05_function_space_geometry.h"
#include "mesh/type_traits.h"

namespace FunctionSpace
{

/** class to get node positions
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceDofsNodesStructured :
  public FunctionSpaceGeometry<MeshType,BasisFunctionType>
{
public:

  //! inherited constructor
  using FunctionSpaceGeometry<MeshType,BasisFunctionType>::FunctionSpaceGeometry;
  
  //! return number of nodes including ghost nodes, i.e. these nodes are known locally but some of them are owned by other ranks
  node_no_t nNodesLocalWithGhosts() const;

  //! return number of nodes in specified coordinate direction
  node_no_t nNodesLocalWithGhosts(int dimension) const;

  //! return number of nodes that are owned by this partition
  node_no_t nNodesLocalWithoutGhosts() const;

  //! return number of nodes in specified coordinate direction that are owned by this partition
  node_no_t nNodesLocalWithoutGhosts(int dimension) const;

  //! return number of dofs
  dof_no_t nDofsLocalWithGhosts() const;
  
  //! return number of dofs
  dof_no_t nDofsLocalWithoutGhosts() const;
  
  //! return number of nodes in specified coordinate direction for the whole global domain
  global_no_t nNodesGlobal(int dimension) const;

  //! return global number of nodes
  global_no_t nNodesGlobal() const;

  //! return global number of dofs
  global_no_t nDofsGlobal() const;

  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;
  
protected:

  // if BasisFunction is Hermite, set the derivative values in the geometry field to mean distance between nodes in the corresponding (x,y,z)-direction
  void setHermiteDerivatives();
};

}  // namespace

#include "function_space/06_function_space_dofs_nodes_structured.tpp"
