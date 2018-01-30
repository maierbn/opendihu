#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.h"

namespace SpatialDiscretization
{

/** general template for any mesh
 */
template<typename BasisOnMeshType, typename IntegratorType, typename Term, typename=typename BasisOnMeshType::Mesh>
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

/** specialisation for Deformable mesh of any dimension D (do proper integration)
 */
template<typename BasisOnMeshType, typename IntegratorType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, IntegratorType, Term, Mesh::isDeformable<typename BasisOnMeshType::Mesh>
> :
  public AssembleStiffnessMatrix<BasisOnMeshType, IntegratorType, Term>
{
public:
  // use constructor of base class
  using AssembleStiffnessMatrix<BasisOnMeshType, IntegratorType, Term>::AssembleStiffnessMatrix;
};

 
};  // namespace

#include "spatial_discretization/finite_element_method/02_stiffness_matrix_stencils.tpp"
