#pragma once

#include "spatial_discretization/finite_element_method/01_matrix.h"

namespace SpatialDiscretization
{

/**
 * Class that prepares the system to enforce Dirichlet boundary conditions, regular numerical integration
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename Dummy= Term>
class BoundaryConditions :
  public FiniteElementMethodMatrix<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodMatrix;

protected:

  //! apply dirichlet boundary conditions, this calls applyBoundaryConditionsWeakForm
  virtual void applyBoundaryConditions();

  //! apply boundary conditions in weak form by adding a term to the rhs
  void applyBoundaryConditionsWeakForm();
};

/**
 * Class that prepares the system to enforce Dirichlet boundary conditions, when Quadrature::None is given, uses values from stiffness matrix
 */
template<typename BasisOnMeshType, typename Term>
class BoundaryConditions<BasisOnMeshType, Quadrature::None, Term, Term> :
  public FiniteElementMethodMatrix<BasisOnMeshType, Quadrature::None, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<BasisOnMeshType, Quadrature::None, Term>::FiniteElementMethodMatrix;

protected:

  //! apply dirichlet boundary conditions, this calls applyBoundaryConditionsWeakForm
  virtual void applyBoundaryConditions();

  //! Apply Dirichlet BC in strong form by setting columns and rows in stiffness matrix to zero such that Dirichlet boundary conditions are met, sets some rows/columns to 0 and the diagonal to 1, changes rhs accordingly
  void applyBoundaryConditionsStrongForm();
};

/**
 * Partial specialization for solid mechanics, mixed formulation
 */
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename QuadratureType,typename Term>
class BoundaryConditions<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, QuadratureType, Term, Equation::isSolidMechanics<Term>> :
  public FiniteElementMethodMatrix<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>, QuadratureType, Term>::FiniteElementMethodMatrix;

protected:
  //! apply dirichlet boundary conditions
  void applyBoundaryConditions(){}

};

};  // namespace

#include "spatial_discretization/finite_element_method/02_boundary_conditions_strong_form.tpp"
#include "spatial_discretization/finite_element_method/02_boundary_conditions_weak_form.tpp"
