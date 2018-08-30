#pragma once

#include "spatial_discretization/finite_element_method/03_assemble_rhs.h"

namespace SpatialDiscretization
{
/** base class implementing right hand side, that can be set by user for poisson equation
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodRhs :
  public AssembleRightHandSide<BasisOnMeshType, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using AssembleRightHandSide<BasisOnMeshType, QuadratureType, Term>::AssembleRightHandSide;

  friend class StiffnessMatrixTester;    ///< a class used for testing
protected:

  //! read in rhs values from config and transfer to weak form
  void setRightHandSide();
};

};  // namespace

#include "spatial_discretization/finite_element_method/04_rhs.tpp"
