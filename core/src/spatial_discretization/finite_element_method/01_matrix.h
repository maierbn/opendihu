#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** Set stiffness and mass matrices with normal integration
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term, typename=typename FunctionSpaceType::Mesh, typename=Term, typename=typename FunctionSpaceType::BasisFunction>
class FiniteElementMethodMatrix :
  public FiniteElementMethodInitializeData<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodInitializeData<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodInitializeData;

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
  public FiniteElementMethodInitializeData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodInitializeData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodInitializeData;

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
  public FiniteElementMethodInitializeData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodInitializeData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodInitializeData;

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
  public FiniteElementMethodInitializeData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodInitializeData<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodInitializeData;

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
  public FiniteElementMethodInitializeData<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodInitializeData<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodInitializeData;

protected:
  //! set entries in mass matrix
  void setStiffnessMatrix();

  //! set entries in mass matrix
  void setMassMatrix();
};*/


/** compute inverse lumped mass matrix
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
class FiniteElementMethodMatrixInverseLumpedMass :
  public FiniteElementMethodMatrix<FunctionSpaceType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodMatrix<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodMatrix;

  //! set entries in the inverse lumped mass matrix, compute from mass matrix
  void setInverseLumpedMassMatrix();
};

extern bool outputAssemble3DStiffnessMatrix_;   ///< if the message about assembly of the 3D stiffness matrix was printed


} // namespace

#include "spatial_discretization/finite_element_method/01_stiffness_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/01_mass_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"
#include "spatial_discretization/finite_element_method/01_stiffness_matrix_integrate.tpp"
#include "spatial_discretization/finite_element_method/01_mass_matrix_integrate.tpp"
#include "spatial_discretization/finite_element_method/01_inverse_lumped_mass_matrix.tpp"
