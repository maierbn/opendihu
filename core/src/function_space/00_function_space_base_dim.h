#pragma once

#include "control/types.h"

namespace FunctionSpace
{

/** base class that provides constants for the numbers of dofs, elements and nodes
  */
template<int D,typename BasisFunctionType>
class FunctionSpaceBaseDim
{
public:

  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerElement();

  //! number of nodes per element
  static constexpr int nNodesPerElement();

  //! number of dofs per node
  static constexpr int nDofsPerNode();

  //! if one assigns every dof to an element it is contained in, the number of degrees of freedom per element (not considering boundary elements)
  static constexpr int averageNDofsPerElement();

  //! if one assigns every node to an element it is contained in, the number of nodes per element (not considering boundary elements)
  static constexpr int averageNNodesPerElement();
};

} // namespace

#include "function_space/00_function_space_base_dim.tpp"
