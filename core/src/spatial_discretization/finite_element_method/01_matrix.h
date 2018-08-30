#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** Set stiffness and mass matrices with normal integration
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term, typename=typename FunctionSpaceType::Mesh, typename=Term, typename=typename FunctionSpaceType::BasisFunction>
class FiniteElementMethodMatrix :
  public FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in stiffness matrix by normal integration
  void setStiffnessMatrix();

  //! set entries in mass matrix by normal integration
  void setMassMatrix();
};

/** stencils
 *  partial specialisation for linear Lagrange, RegularFixed mesh, dimension 1 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<1ul>,
  Equation::hasLaplaceOperator<Term>,
  BasisFunction::LagrangeOfOrder<1>
> :
  public FiniteElementMethodBase<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in stiffness matrix by using stencils
  void setStiffnessMatrix();

  //! set entries in mass matrix by using stencils
  void setMassMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 2 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<2ul>,
  Equation::hasLaplaceOperator<Term>,
  BasisFunction::LagrangeOfOrder<1>
> :
  public FiniteElementMethodBase<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in stiffness matrix by using stencils
  void setStiffnessMatrix();

  //! set entries in mass matrix by using stencils
  void setMassMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 3 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<3ul>,
  Equation::hasLaplaceOperator<Term>,
  BasisFunction::LagrangeOfOrder<1>
> :
  public FiniteElementMethodBase<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in stiffness matrix by using stencils
  void setStiffnessMatrix();

  //! set entries in mass matrix by using stencils
  void setMassMatrix();
};

/** specialisation for Deformable mesh of any dimension D (do proper integration)
 *
 * Class that creates the stiffness matrix and mass matrix by integrating its integrand over the elements.
 * What to integrate is given by the class template Term.
 */
/*template<typename FunctionSpaceType, typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  FunctionSpaceType,
  QuadratureType,
  Term,
  Equation::doesNotUseStencilsNorSolidMechanics<typename FunctionSpaceType::BasisFunction,typename FunctionSpaceType::Mesh,Term>,
  Term,
  typename FunctionSpaceType::BasisFunction
> :
  public FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in mass matrix
  void setStiffnessMatrix();

  //! set entries in mass matrix
  void setMassMatrix();
};*/



};  // namespace

#include "spatial_discretization/finite_element_method/01_stiffness_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/01_mass_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"
#include "spatial_discretization/finite_element_method/01_stiffness_matrix_integrate.tpp"
#include "spatial_discretization/finite_element_method/01_mass_matrix_integrate.tpp"
