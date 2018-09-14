#pragma once

#include <memory>

#include "spatial_discretization/finite_element_method/01_matrix.h"
#include "spatial_discretization/dirichlet_boundary_conditions.h"

namespace SpatialDiscretization
{

/**
 * Class that prepares the system to enforce Dirichlet boundary conditions.
 * It uses values from the already set up system matrix and transfers them to rhs.
 * In this way it works for stencil based matrices as well as matrices created by regular numerical integration.
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term, typename Dummy= Term>
class BoundaryConditions :
  public FiniteElementMethodMatrix<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodMatrix;

protected:

  //! apply dirichlet boundary conditions, this calls applyBoundaryConditionsWeakForm
  virtual void applyBoundaryConditions();

  //! parse config and fill local member variables
  void parseBoundaryConditions();

  std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions_ = nullptr;  ///< object that parses Dirichlet boundary conditions and applies them to system matrix and rhs
};

/**
 * Partial specialization for solid mechanics, mixed formulation
 */
template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename QuadratureType,typename Term>
class BoundaryConditions<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, Term, Equation::isSolidMechanics<Term>> :
  public FiniteElementMethodMatrix<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, Term>::FiniteElementMethodMatrix;

protected:
  //! apply dirichlet boundary conditions
  void applyBoundaryConditions(){}

};

};  // namespace

#include "spatial_discretization/finite_element_method/02_boundary_conditions.tpp"
