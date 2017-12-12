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
    
  //! fill a vector with node positions 
  void getNodePositions(std::vector<double> &nodePositions);
  
  friend class NodePositionsTester;
private:
 
  Vec nodePositions_;                         ///< contains the positions of the nodes in physical space, for 2D consecutive (x,y) pairs, for 3D consecutive (x,y,z) triples
};  

}  // namespace

#include "mesh/deformable.tpp"