#pragma once

#include <array>
#include "control/types.h"

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
class BasisOnMeshBaseDim
{
public:
 
  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerElement();
  
  //! number of nodes per element
  static constexpr int nNodesPerElement();
  
  //! number of dofs per node
  static constexpr int nDofsPerNode();
  
  //! if one assigns every dof to an element it is contained in, the number of degrees of freedom per element (not considering border elements)
  static constexpr int averageNDofsPerElement();
  
  //! if one assigns every node to an element it is contained in, the number of nodes per element (not considering border elements)
  static constexpr int averageNNodesPerElement();
};
 
};  // namespace

#include "basis_on_mesh/00_basis_on_mesh_base_dim.tpp"