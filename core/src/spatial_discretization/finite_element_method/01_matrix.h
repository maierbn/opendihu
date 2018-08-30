#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** Set stiffness and mass matrices with normal integration
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename=typename BasisOnMeshType::Mesh, typename=Term, typename=typename BasisOnMeshType::BasisFunction>
class FiniteElementMethodMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;

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
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<1ul>,
  Equation::hasLaplaceOperator<Term>,
  BasisFunction::LagrangeOfOrder<1>
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;

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
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<2ul>,
  Equation::hasLaplaceOperator<Term>,
  BasisFunction::LagrangeOfOrder<1>
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;

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
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<3ul>,
  Equation::hasLaplaceOperator<Term>,
  BasisFunction::LagrangeOfOrder<1>
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;

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
/*template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Equation::doesNotUseStencilsNorSolidMechanics<typename BasisOnMeshType::BasisFunction,typename BasisOnMeshType::Mesh,Term>,
  Term,
  typename BasisOnMeshType::BasisFunction
> :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;

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
