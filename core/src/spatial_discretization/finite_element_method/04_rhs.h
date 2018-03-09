#pragma once

#include "spatial_discretization/finite_element_method/02_stiffness_matrix.h"
#include "spatial_discretization/finite_element_method/03_assemble_rhs.h"

namespace SpatialDiscretization
{
 
/** base class for baseRhs
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename=typename BasisOnMeshType::Mesh, typename=Term>
class FiniteElementMethodBaseRhs :
  public FiniteElementMethodStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>
{
};

/** specialization for linear Lagrange, 1D regular mesh (uses stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<1ul>, Equation::hasLaplaceOperator<Term>> :
  public FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
    ::FiniteElementMethodStiffnessMatrix;
    
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
  
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void setMassMatrix();
};

/** specialization for linear Lagrange, 2D regular mesh (uses stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<2ul>, Equation::hasLaplaceOperator<Term>> :
  public FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
    ::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
  
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void setMassMatrix();
};

/** specialization for linear Lagrange, 3D regular mesh (uses stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodBaseRhs<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term, Mesh::StructuredRegularFixedOfDimension<3ul>, Equation::hasLaplaceOperator<Term>> :
  public FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodStiffnessMatrix<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term>
    ::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
  
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void setMassMatrix();
};

/** specialisation for RegularFixed not linear Lagrange, and other meshes of any dimension D (do proper integration of rhs)
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodBaseRhs<BasisOnMeshType, QuadratureType, Term, 
    Equation::doesNotUseStencils<typename BasisOnMeshType::BasisFunction,typename BasisOnMeshType::Mesh,Term>> :
  public AssembleRightHandSide<BasisOnMeshType, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using AssembleRightHandSide<BasisOnMeshType, QuadratureType, Term>::AssembleRightHandSide;
  
};

/** base class implementing right hand side, that can be set by user for poisson equation
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodRhs :
  public FiniteElementMethodBaseRhs<BasisOnMeshType, QuadratureType, Term>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBaseRhs<BasisOnMeshType, QuadratureType, Term>
    ::FiniteElementMethodBaseRhs;
 
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  //! read in rhs values from config and transfer to weak form
  void setRightHandSide();
};
 
};  // namespace

#include "spatial_discretization/finite_element_method/04_rhs.tpp"
#include "spatial_discretization/finite_element_method/04_rhs_stencils.tpp"
