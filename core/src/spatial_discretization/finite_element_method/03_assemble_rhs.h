#pragma once

#include "spatial_discretization/finite_element_method/02_boundary_conditions.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** specialisation for RegularFixed not linear Lagrange, and other meshes of any dimension D (do proper integration of rhs)
 *
 * Class that sets the right hand side vector by integrating the integrand over the elements.
 * What to integrate is given by the class template Term.
 */
template<typename FunctionSpaceType, typename QuadratureType, int nComponents, typename Term, typename=Term>
class AssembleRightHandSide :
  public BoundaryConditions<FunctionSpaceType, QuadratureType, nComponents, Term>
{
public:
  // use constructor of base class
  using BoundaryConditions<FunctionSpaceType, QuadratureType, nComponents, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRightHandSideWithMassMatrix();
};

/** specialization for linear Lagrange, 1D regular mesh (uses stencils)
 */
template<typename QuadratureType, typename Term>
class AssembleRightHandSide<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term, Equation::hasLaplaceOperator<Term>> :
  public BoundaryConditions<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term>
{
public:
  //! use constructor of base class
  using BoundaryConditions<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRightHandSideWithMassMatrix();
};

/** specialization for linear Lagrange, 2D regular mesh (uses stencils)
 */
template<typename QuadratureType, typename Term>
class AssembleRightHandSide<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term, Equation::hasLaplaceOperator<Term>> :
  public BoundaryConditions<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term>
{
public:
  //! use constructor of base class
  using BoundaryConditions<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRightHandSideWithMassMatrix();
};

/** specialization for linear Lagrange, 3D regular mesh (uses stencils)
 */
template<typename QuadratureType, typename Term>
class AssembleRightHandSide<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term, Equation::hasLaplaceOperator<Term>> :
  public BoundaryConditions<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term>
{
public:
  //! use constructor of base class
  using BoundaryConditions<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, 1, Term>
    ::BoundaryConditions;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRightHandSideWithMassMatrix();
};

/** specialisation for RegularFixed not linear Lagrange, and other meshes of any dimension D (do proper integration of rhs)
 *
 * Class that sets the right hand side vector by integrating the integrand over the elements.
 * What to integrate is given by the class template Term.
 *//*
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
class AssembleRightHandSide<FunctionSpaceType, QuadratureType, 1, Term,
    Equation::doesNotUseStencils<typename FunctionSpaceType::BasisFunction,typename FunctionSpaceType::Mesh,Term>> :
  public BoundaryConditions<FunctionSpaceType, QuadratureType, 1, Term>
{
public:
  // use constructor of base class
  using BoundaryConditions<FunctionSpaceType, QuadratureType, 1, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRightHandSideWithMassMatrix();
};*/

/**
 * Partial specialization for solid mechanics, mixed formulation
 */
/*template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
class AssembleRightHandSide<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, nComponents, Term,
    Equation::isSolidMechanics<Term>> :
  public BoundaryConditions<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, nComponents, Term>
{
public:
  // use constructor of base class
  using BoundaryConditions<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, nComponents, Term>::BoundaryConditions;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRightHandSideWithMassMatrix(){}
};*/

} // namespace

#include "spatial_discretization/finite_element_method/03_assemble_rhs_integrate.tpp"
#include "spatial_discretization/finite_element_method/03_assemble_rhs_stencils.tpp"
