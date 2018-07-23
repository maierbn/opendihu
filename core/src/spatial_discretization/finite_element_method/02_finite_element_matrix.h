#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/01_assemble_finite_element_matrix.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** general template for any mesh, this should not be in use
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename=typename BasisOnMeshType::Mesh, typename=Term, typename=typename BasisOnMeshType::BasisFunction>
class FiniteElementMethodMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix(){LOG(FATAL)<<"FiniteElementMethodStiffnessMatrix inheritance is wrong!";}
  //! set entries in mass matrix
  void setMassMatrix(){LOG(FATAL)<<"FiniteElementMethodMassMatrix inheritance is wrong!";}
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
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
  //! set entries in mass matrix
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
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
  //! set entries in mass matrix
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
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
  //! set entries in mass matrix
  void setMassMatrix();
};

/** specialisation for Deformable mesh of any dimension D (do proper integration)
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Equation::doesNotUseStencilsNorSolidMechanics<typename BasisOnMeshType::BasisFunction,typename BasisOnMeshType::Mesh,Term>,
  Term,
  typename BasisOnMeshType::BasisFunction
> :
  public AssembleFiniteElementMatrix<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using AssembleFiniteElementMatrix<BasisOnMeshType, QuadratureType, Term>::AssembleFiniteElementMatrix;
};




};  // namespace

#include "spatial_discretization/finite_element_method/02_stiffness_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/02_mass_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"
