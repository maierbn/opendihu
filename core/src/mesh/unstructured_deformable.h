#pragma once

#include <vector>
#include <petscmat.h>

#include "control/types.h"
#include "mesh/mesh.h"
#include "mesh/deformable.h"
#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"

namespace Mesh
{

/**
 * An arbitrary mesh where each element can be adjacent to any other element or at the border of the computational domain.
 * There is no restriction that the total domain must be quadratic or cubic.
 */
template<int D>
class UnstructuredDeformable : public MeshD<D>, public Deformable
{
public:
  //! constructor of base class
  using MeshD<D>::MeshD;
  
 
private:
};  

}  // namespace

#include "mesh/unstructured_deformable.tpp"