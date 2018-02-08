#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.h"
#include "equation/solid_mechanics.h"

namespace SpatialDiscretization
{

/** general template for any mesh
 */
template<typename BasisOnMeshType, typename IntegratorType, typename Term, typename=typename BasisOnMeshType::Mesh, typename=Term>
class FiniteElementMethodStiffnessMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, IntegratorType>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

/** stencils
 *  partial specialisation for linear Lagrange, RegularFixed mesh, dimension 1 (use precomputed stencils)
 */
template<typename IntegratorType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>>,IntegratorType, Term
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>>, IntegratorType>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>>, IntegratorType>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 2 (use precomputed stencils)
 */
template<typename IntegratorType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>>, IntegratorType, Term
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>>, IntegratorType>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>>, IntegratorType>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 3 (use precomputed stencils)
 */
template<typename IntegratorType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>>, IntegratorType, Term
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>>, IntegratorType>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>>, IntegratorType>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};
 
/** specialisation for solid mechanics
 */
template<typename BasisOnMeshType, typename MixedIntegratorType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, MixedIntegratorType, Term, Mesh::isDeformable<typename BasisOnMeshType::Mesh>, Equation::isSolidMechanics<Term>
> :
  public FiniteElementMethodBase<BasisOnMeshType, MixedIntegratorType>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, MixedIntegratorType>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! no set right hand side functionality
  void setRightHandSide(){}

private:  
  //! compute the deformation gradient from geometry and jacobian data of an element at parameter xi
  std::array<Vec3,BasisOnMeshType::dim()> computeDeformationGradient(std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &geometryReference,
                                                                     std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &displacement, 
                                                                     std::array<Vec3,BasisOnMeshType::dim()> &jacobian, 
                                                                     std::array<double, BasisOnMeshType::dim()> xi);
};

/** specialisation for Deformable mesh of any dimension D (do proper integration)
 */
template<typename BasisOnMeshType, typename IntegratorType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, IntegratorType, Term, Mesh::isDeformable<typename BasisOnMeshType::Mesh>, Equation::hasLaplaceOperator<Term>
> :
  public AssembleStiffnessMatrix<BasisOnMeshType, IntegratorType, Term>
{
public:
  // use constructor of base class
  using AssembleStiffnessMatrix<BasisOnMeshType, IntegratorType, Term>::AssembleStiffnessMatrix;  
};

 
};  // namespace

#include "spatial_discretization/finite_element_method/02_stiffness_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/02_stiffness_matrix_solid_mechanics.tpp"
