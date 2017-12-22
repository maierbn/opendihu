#pragma once

#include <array>

#include "mesh/mesh.h"
#include "control/types.h"

#include <petscmat.h>

#include "mesh/regular.h"

namespace Mesh
{

template<unsigned long D>
class Deformable : public Regular<D>
{
public:
 
  //! construct mesh from python settings
  Deformable(PyObject *specificSettings);
    
  //! fill a vector with geometry field entries (node position for Lagrange elements) 
  void getGeometry(std::vector<double> &geometry);
  
  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_idx_t dofNo);
  
  friend class NodePositionsTester;
private:
 
  Vec nodePositions_;  ///< contains the positions of the nodes in physical space, always as consecutive (x,y,z) triple. For 1D/2D only problems, y,z resp. z are set to 0.
};  

}  // namespace

#include "mesh/deformable.tpp"