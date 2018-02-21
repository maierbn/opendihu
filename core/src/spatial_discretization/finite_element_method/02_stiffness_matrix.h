#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.h"
#include "equation/solid_mechanics.h"
#include "equation/compressible_mooney_rivlin.h"

namespace SpatialDiscretization
{

/** general template for any mesh
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename=typename BasisOnMeshType::Mesh, typename=Term>
class FiniteElementMethodStiffnessMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

/** stencils
 *  partial specialisation for linear Lagrange, RegularFixed mesh, dimension 1 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>,QuadratureType, Term
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 2 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<2ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

/** partial specialisation for linear Lagrange, RegularFixed mesh, dimension 3 (use precomputed stencils)
 */
template<typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType, Term
> :
  public FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType>
{
public:
  //! use constructor of base class
  using FiniteElementMethodBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<3ul>, BasisFunction::LagrangeOfOrder<1>>, QuadratureType>::FiniteElementMethodBase;
  
protected:
  //! set entries in stiffness matrix
  void setStiffnessMatrix();
};

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
class FiniteElementMethodSolidMechanicsUtility
{
public:

protected:  
  //! compute the deformation gradient from geometry and jacobian data of an element at parameter xi
  std::array<Vec3,BasisOnMeshType::dim()> computeDeformationGradient(const std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &geometryReference,
                                                                     const std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &displacement, 
                                                                     const std::array<Vec3,BasisOnMeshType::dim()> &jacobian, 
                                                                     const std::array<double, BasisOnMeshType::dim()> xi);

  //! compute right cauchy green tensor
  std::array<Vec3,BasisOnMeshType::dim()> computeRightCauchyGreenTensor(const std::array<Vec3,BasisOnMeshType::dim()> &deformationGradient);

  //! compute the standard invariants I1=tr(C), I2=1/2(tr(C)^2 - tr(C^2)), I3=det(C)
  std::array<double,3> computeInvariants(const std::array<Vec3,BasisOnMeshType::dim()> &rightCauchyGreen, const double determinant);

  
  //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC
  std::array<Vec3,3> computePK2Stress(const std::array<Vec3,3> &rightCauchyGreen,
                                      const std::array<Vec3,3> &inverseRightCauchyGreen,
                                      const std::array<double,3> invariants);
  
  //! compute the elasticity tensor C = 2*dS/dC. Due to hyperelasticity there are symmetries C_{ijrs} = C_{jirs} and C_{ijrs} = C_{rsij} that leave 21 independent values.
  std::array<double,21> computeElasticityTensor(const std::array<Vec3,3> &rightCauchyGreen,
                                                const std::array<Vec3,3> &inverseRightCauchyGreen,
                                                const std::array<double,3> invariants);
  
  //! compute the pressure from displacements, using the formula 2.8 of Sussman and Bathe "A finite element formulation for nonlinear incompressible elastic and inelastic analysis"
  double computePressureFromDisplacements(double deformationGradientDeterminant, const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &PK2Stress);
};
 
/** specialisation for incompressible solid mechanics
 */
template<typename BasisOnMeshType, typename MixedQuadratureType>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, MixedQuadratureType, Equation::Static::SolidMechanics, Mesh::isDeformable<typename BasisOnMeshType::Mesh>
> :
  public FiniteElementMethodBase<BasisOnMeshType, MixedQuadratureType>,
  public FiniteElementMethodSolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Equation::Static::SolidMechanics>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, MixedQuadratureType>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! no set right hand side functionality
  void setRightHandSide(){}
};
 
/** specialisation for compressible mooney rivlin
 */
template<typename BasisOnMeshType, typename MixedQuadratureType>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, MixedQuadratureType, Equation::Static::CompressibleMooneyRivlin, Mesh::isDeformable<typename BasisOnMeshType::Mesh>
> :
  public FiniteElementMethodBase<BasisOnMeshType, MixedQuadratureType>,
  public FiniteElementMethodSolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Equation::Static::CompressibleMooneyRivlin>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, MixedQuadratureType>::FiniteElementMethodBase;
  
  //! assemble the stiffness matrix
  void setStiffnessMatrix();
  
  //! no set right hand side functionality
  void setRightHandSide(){}
  
  //! initialize material constants from config
  void initialize() override;
  
private:
  
  //! compute the reduced invariants J1 = I1*I3^-1/3, J2 = I2*I3^-2/3, J3=det F
  std::array<double,3> computeReducedInvariants(const std::array<double,3> &invariants, double deformationGradientDeterminant);

  //! compute the elasticity tensor C = 2*sym(dS/dC). Due to hyperelasticity there are symmetries C_{ijrs} = C_{jirs} and C_{ijrs} = C_{rsij} that leave 21 independent values.
  std::array<double,21> computeElasticityTensor(const std::array<Vec3,3> &rightCauchyGreen, 
                                                const std::array<Vec3,3> &inverseRightCauchyGreen, 
                                                const std::array<double,3> invariants, 
                                                const std::array<double,3> &reducedInvariants);
  
  //! compute 2nd Piola-Kirchhoff stress tensor S = 2*sym(dPsi/dC)
  std::array<Vec3,3> computePK2Stress(const std::array<Vec3,3> &rightCauchyGreen, 
                                      const std::array<Vec3,3> &inverseRightCauchyGreen, 
                                      const std::array<double,3> invariants, 
                                      const std::array<double,3> &reducedInvariants);
 
  //! return the index to the elastcity array with 21 distinct entries, for logical component klrs
  int getElasticityEntryNo(int k, int l, int r, int s);
  
  //! return the entry klrs of the elasticity tensor
  double getElasticityEntry(std::array<double, 21> &elasticity, int k, int l, int r, int s);
  
  double c1_, c2_;  ///< material constants for mooney-rivlin material
  double kappa_;  ///< bulk modulus
 
};

/** specialisation for Deformable mesh of any dimension D (do proper integration)
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, QuadratureType, Term, Mesh::isDeformable<typename BasisOnMeshType::Mesh>, Equation::hasLaplaceOperator<Term>
> :
  public AssembleStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using AssembleStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>::AssembleStiffnessMatrix;  
};

 
};  // namespace

#include "spatial_discretization/finite_element_method/02_stiffness_matrix_stencils.tpp"
#include "spatial_discretization/finite_element_method/02_stiffness_matrix_solid_mechanics.tpp"
