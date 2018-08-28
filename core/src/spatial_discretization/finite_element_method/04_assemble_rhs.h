#pragma once

#include "spatial_discretization/finite_element_method/03_boundary_conditions.h"

namespace SpatialDiscretization
{

/**
 * Class that sets the right hand side vector by integrating the integrand over the elements.
 * What to integrate is given by the class template Term.
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename Dummy=Term>
class AssembleRightHandSide :
  public BoundaryConditions<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using BoundaryConditions<BasisOnMeshType, QuadratureType, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();

  //! set the matrix that transforms a vector of rhs value into a vector that contains the right hand side in discretized form for FEM
  void setMassMatrix();
};

/**
 * Partial specialization for solid mechanics, mixed formulation
 */
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename QuadratureType,typename Term>
class AssembleRightHandSide<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, QuadratureType, Term, Equation::isSolidMechanics<Term>> :
public BoundaryConditions<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, QuadratureType, Term>
{
public:
  // use constructor of base class
  using BoundaryConditions<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, QuadratureType, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm(){}

  //! set the matrix that transforms a vector of rhs value into a vector that contains the right hand side in discretized form for FEM
  void setMassMatrix(){}
};

};  // namespace

#include "spatial_discretization/finite_element_method/04_assemble_rhs.tpp"
