#pragma once

#include <memory>

#include "spatial_discretization/finite_element_method/01_matrix.h"
#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"
#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions.h"

namespace SpatialDiscretization
{

/**
 * Class that prepares the system to enforce Dirichlet boundary conditions.
 * It uses values from the already set up system matrix and transfers them to rhs.
 * In this way it works for stencil based matrices as well as matrices created by regular numerical integration.
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term, typename Dummy= Term>
class BoundaryConditions :
  public FiniteElementMethodMatrixInverseLumpedMass<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrixInverseLumpedMass<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodMatrixInverseLumpedMass;

  //! enable or disable boundary condition handling on initialization, set to false to not care for boundary conditions
  virtual void setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled);

  //! set the dirichlet boundary condition object, this can be done to set BC from within the C++ code, not from python config
  void setDirichletBoundaryConditions(std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions);

  //! set the neumann boundary condition object, this can be done to set BC from within the C++ code, not from python config
  void setNeumannBoundaryConditions(std::shared_ptr<NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,1>> neumannBoundaryConditions);

protected:

  //! apply dirichlet boundary conditions, this calls applyBoundaryConditionsWeakForm
  virtual void applyBoundaryConditions();

  //! parse config and fill local member variables
  void parseBoundaryConditions();

  bool boundaryConditionHandlingEnabled_ = true;   ///< if the boundary conditions should be handled in this class, if false, nothing is done here. This is the case if the FiniteElementMethod is used within a timestepping scheme. Then the time stepping scheme constructs its system matrix out of this class' stiffness matrix and applied Dirichlet boundary condition handle there.
  std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions_ = nullptr;  ///< object that parses Dirichlet boundary conditions and applies them to system matrix and rhs
  std::shared_ptr<NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,1>> neumannBoundaryConditions_ = nullptr;  ///< object that parses Neumann boundary conditions and applies them to the rhs
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

} // namespace

#include "spatial_discretization/finite_element_method/02_boundary_conditions.tpp"
