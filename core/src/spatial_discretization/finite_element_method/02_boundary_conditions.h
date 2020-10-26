#pragma once

#include <memory>

#include "spatial_discretization/finite_element_method/01_matrix.h"
#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.h"
#include "spatial_discretization/neumann_boundary_conditions/01_neumann_boundary_conditions.h"

namespace SpatialDiscretization
{

/**
 * Class that prepares the system to enforce Dirichlet boundary conditions.
 * It uses values from the already set up system matrix and transfers them to rhs.
 * In this way it works for stencil based matrices as well as matrices created by regular numerical integration.
 */
template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term,typename Dummy= Term>
class BoundaryConditions :
  public FiniteElementMethodMatrixInverseLumpedMass<FunctionSpaceType,QuadratureType,nComponents,Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrixInverseLumpedMass<FunctionSpaceType,QuadratureType,nComponents,Term>::FiniteElementMethodMatrixInverseLumpedMass;

  //! enable or disable boundary condition handling on initialization, set to false to not care for boundary conditions
  virtual void setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled);

  //! set the dirichlet boundary condition object, this can be done to set BC from within the C++ code, not from python config
  void setDirichletBoundaryConditions(std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,nComponents>> dirichletBoundaryConditions);

  //! set the neumann boundary condition object, this can be done to set BC from within the C++ code, not from python config
  void setNeumannBoundaryConditions(std::shared_ptr<NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents>> neumannBoundaryConditions);

  //! reset to pre-initialized state
  virtual void reset();

  //! apply dirichlet boundary conditions, this calls applyBoundaryConditionsWeakForm
  virtual void applyBoundaryConditions();

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> debuggingFieldVariable_;  //< temporary debugging field variable that will be output
protected:

  //! parse config and fill local member variables
  void parseBoundaryConditions();

  //! apply the dirichlet type boundary conditions
  void applyDirichletBoundaryConditions();

  //! apply the neumann type boundary conditions
  void applyNeumannBoundaryConditions();

  bool boundaryConditionHandlingEnabled_ = true;    //< if the boundary conditions should be handled in this class, if false, nothing is done here. This is the case if the FiniteElementMethod is used within a timestepping scheme. Then the time stepping scheme constructs its system matrix out of this class' stiffness matrix and applied Dirichlet boundary condition handle there.
  std::shared_ptr<DirichletBoundaryConditions<FunctionSpaceType,nComponents>> dirichletBoundaryConditions_ = nullptr;             //< object that parses Dirichlet boundary conditions and applies them to system matrix and rhs
  std::shared_ptr<NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents>> neumannBoundaryConditions_ = nullptr;  //< object that parses Neumann boundary conditions and applies them to the rhs
  bool systemMatrixAlreadySet_ = false;             //< if the system matrix has been changed to account for dirichlet DCs, which means that rows/columns of BC dofs were set to zero and diagonal to 1
  bool neumannBoundaryConditionsApplied_ = false;   //< if the neumann BC were already applied
  bool dirichletBoundaryConditionsApplied_ = false; //< if the dirichlet BC were already applied after the last initialize() or setDirichletBoundaryConditions()

};

/**
 * Partial specialization for solid mechanics, mixed formulation
 */
/*template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
class BoundaryConditions<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>, QuadratureType, nComponents, Term, Equation::isSolidMechanics<Term>> :
  public FiniteElementMethodMatrix<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,QuadratureType,nComponents,Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,QuadratureType,nComponents,Term>::FiniteElementMethodMatrix;

protected:
  //! apply dirichlet boundary conditions
  void applyBoundaryConditions(){}

};*/

} // namespace

#include "spatial_discretization/finite_element_method/02_boundary_conditions.tpp"
