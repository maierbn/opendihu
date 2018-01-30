#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "control/types.h"
#include "mesh/mesh.h"

namespace Mesh
{
  
/**
 * A deformable mesh (interface), i.e. it can change its geometry while computation. This is only an interface.
 */
class Deformable
{
public:
 
  //! fill a vector with positions of the nodes, consecutive (x,y,z) values
  virtual void getNodePositions(std::vector<double> &nodePositions) const = 0;
  
  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  virtual Vec3 getGeometry(node_idx_t dofNo) const = 0;
};

};    // namespace
