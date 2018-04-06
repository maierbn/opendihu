#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.h"
#include "equation/type_traits.h"

namespace SpatialDiscretization
{

/** general template for any mesh
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename=typename BasisOnMeshType::Mesh, typename=Term>
class FiniteElementMethodStiffnessMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** stencils
 *  partial specialisation for linear Lagrange, RegularFixed mesh, dimension 1 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType, 
  Term, 
  Mesh::StructuredRegularFixedOfDimension<1ul>, 
  Equation::hasLaplaceOperator<Term>
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 2 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, 
  QuadratureType, 
  Term, 
  Mesh::StructuredRegularFixedOfDimension<2ul>,
  Equation::hasLaplaceOperator<Term>
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 3 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>,
  QuadratureType,
  Term,
  Mesh::StructuredRegularFixedOfDimension<3ul>,
  Equation::hasLaplaceOperator<Term>
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** specialisation for Deformable mesh of any dimension D (do proper integration)
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, QuadratureType, Term, Equation::doesNotUseStencilsNorSolidMechanics<typename BasisOnMeshType::BasisFunction,typename BasisOnMeshType::Mesh, Term>
> :
  public AssembleStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using AssembleStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>::AssembleStiffnessMatrix;  
};

 
};  // namespace

#include "spatial_discretization/finite_element_method/02_stiffness_matrix_stencils.tpp"
