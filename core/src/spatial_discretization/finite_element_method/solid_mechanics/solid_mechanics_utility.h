#pragma once

#include <Python.h>  // has to be the first included header

namespace SpatialDiscretization
{
 
/** helper class that encapsulates functionality for solid mechanics
 */
template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
class SolidMechanicsUtility
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
 
};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.tpp"
