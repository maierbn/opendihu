#pragma once

#include <Python.h>  // has to be the first included header

#include <petscmat.h>
#include <iostream>
#include <memory>
#include <map>

#include "control/types.h"
#include "field_variable/exfile_representation.h"
#include "field_variable/element_to_node_mapping.h"
#include "field_variable/node_to_dof_mapping.h"

namespace FieldVariable
{

/** For every element the adjacent dofs.
 *  Also for global nodes the dofs are stored, this is needed for the construction of the dof mapping.
 */
class ElementToDofMapping
{
public:
 
  //! resize internal representation variable to number of elements
  void setNumberElements(element_idx_t nElements);
 
  //! setup the element to dof mapping and create a node to dof mapping on the fly
  std::shared_ptr<NodeToDofMapping> setup(std::shared_ptr<ExfileRepresentation> exfileRepresentation,
                                          std::shared_ptr<ElementToNodeMapping> elementToNodeMapping,
                                          const int nDofsPerNode);
  
  //! get all dofs of an element
  std::vector<int> &getElementDofs(element_idx_t elementGlobalNo);
  
  //! return the number of dofs
  int nDofs() const;
  
  //! comparison operator
  bool operator==(const ElementToDofMapping &rhs);
  
private:
  std::vector<std::vector<int>> dofs_;  ///< for every element the list of dofs
  int nDofs_ = 0;            ///< total number of dofs
};

};  // namespace
