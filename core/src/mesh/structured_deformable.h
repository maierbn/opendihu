#pragma once

#include <vector>
#include <petscmat.h>

#include "control/types.h"
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "mesh/deformable.h"

namespace Mesh
{

/**
 * A structured mesh, i.e. fixed number of elements in each coordinate direction. 
 * The shape of the elements/the node positions is arbitrary.
 * The node positions can be changed by computation, e.g. for computational mechanics.
 */
template<int D>
class StructuredDeformable : public Structured<D>, public Deformable
{
public: 
  //! constructor of base class
  using Structured<D>::Structured;
  
};  

}  // namespace

#include "mesh/structured_deformable.tpp"