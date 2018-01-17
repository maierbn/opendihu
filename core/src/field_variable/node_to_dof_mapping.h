#pragma once

#include <petscmat.h>
#include <iostream>
#include <memory>
#include <map>

#include "control/types.h"
#include "field_variable/exfile_representation.h"
#include "field_variable/element_to_node_mapping.h"

namespace FieldVariable
{

/** For every element the adjacent dofs.
 *  Also for global nodes the dofs are stored, this is needed for the construction of the dof mapping.
 */
class NodeToDofMapping
{
public:
 
  /** all the dofs of a node, also the indices at which position of the node in the exfiles the value occurs
   */
  struct NodeDofInformation
  {
    std::vector<int> dofs;               ///< the dofs of the node
    std::vector<double> scaleFactors;    ///< the scale factors for Hermite basis functions, not used yet, but set maybe for later use
    std::vector<int> exfileValueIndices; ///< the value indices in the exnode file, entry i of dofs[i] and exfileValueIndices[i] correspond
  };
  
  //! get the dof information of a node
  NodeDofInformation &getNodeDofInformation(node_idx_t nodeGlobalNo);
  
  //! get the dofs of a node
  std::vector<int> &getNodeDofs(node_idx_t nodeGlobalNo);
  
  //! get the dofs of a node
  std::vector<double> &getNodeScaleFactors(node_idx_t nodeGlobalNo);
  
  //! return the number of nodes
  int nNodes() const;
  
  //! check if a node is already contained in the internal node dof information map
  bool containsNode(node_idx_t nodeGlobalNo) const;
  
  //! return the number of versions of a particular node. A version is an OpenCMISS iron construct that allows multiple dofs of the same type on a single node, useful for modeling discontinuities
  int getNumberVersions(node_idx_t nodeGlobalNo, const int nNodesPerElement);
  
private:
  std::map<node_idx_t, NodeDofInformation> nodeDofInformation_;   ///< for global node number the associated dofs, node-to-dof mapping
};

};  // namespace