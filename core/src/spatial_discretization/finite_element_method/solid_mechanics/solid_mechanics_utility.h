#pragma once

#include <Python.h>  // has to be the first included header

#include "spatial_discretization/finite_element_method/solid_mechanics/elasticity_tensor.h"

namespace SpatialDiscretization
{
 
/** helper class that encapsulates functionality for solid mechanics
 */
template<typename BasisOnMeshType, typename Term>
class SolidMechanicsUtility
{
public:

protected:
  //! compute the deformation gradient, F = du/dX, from the displacements and jacobian data of an element at point in parameter space xi. Note that the derivative is w.r.t X, i.e. reference configuration and not parameter space.
  std::array<Vec3,BasisOnMeshType::dim()> computeDeformationGradient(const std::array<Vec3,BasisOnMeshType::nDofsPerElement()> &displacement,
                                                                     const std::array<Vec3,BasisOnMeshType::dim()> &inverseJacobianMaterial,
                                                                     const std::array<double, BasisOnMeshType::dim()> xi);
  
  //! compute the deformation gradient w.r.t parameter space, F = du/dxi, at point xi in parameter space. This is not used for FE, because there we need the deformation gradient w.r.t X
  std::array<Vec3,BasisOnMeshType::dim()> computeDeformationGradientParameterSpace(const std::array<Vec3,BasisOnMeshType::nDofsPerElement()> &displacement,
                                                                                   const std::array<double, BasisOnMeshType::dim()> xi); 
  //! compute right cauchy green tensor, C
  std::array<Vec3,BasisOnMeshType::dim()> computeRightCauchyGreenTensor(const std::array<Vec3,BasisOnMeshType::dim()> &deformationGradient);

  //! compute the standard invariants I1=tr(C), I2=1/2(tr(C)^2 - tr(C^2)), I3=det(C)
  std::array<double,3> computeInvariants(const std::array<Vec3,BasisOnMeshType::dim()> &rightCauchyGreen, 
                                         const double rightCauchyGreenDeterminant);

  //! compute the reduced invariants Ibar1=J^{-2/3}*I1, Ibar2=J^{-4/3}*I2
  std::array<double,2> computeReducedInvariants(const std::array<double,3> invariants, 
                                                const double deformationGradientDeterminant);
  
  //! compute the artifical pressure for penalty formulation, p = dPsi_vol/dJ, pTilde = p + J*dp/dJ
  double computeArtificialPressure(const double deformationGradientDeterminant, 
                                   double &artificialPressureTilde);
  
  //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
  std::array<Vec3,3> computePK2Stress(const double pressure, 
                                      const std::array<Vec3,3> &rightCauchyGreen,
                                      const std::array<Vec3,3> &inverseRightCauchyGreen,
                                      std::array<double,2> reducedInvariants,
                                      const double deformationGradientDeterminant,
                                      std::array<Vec3,3> &fictitiousPK2Stress,
                                      std::array<Vec3,3> &pk2StressIsochoric
                                     );
  
  //! compute the elasticity tensor C = 2*dS/dC. Due to hyperelasticity there are symmetries C_{ijrs} = C_{jirs} and C_{ijrs} = C_{rsij} that leave 21 independent values.
  ElasticityTensor computeElasticityTensorCoupledStrainEnergy(const std::array<Vec3,3> &rightCauchyGreen,
                                                              const std::array<Vec3,3> &inverseRightCauchyGreen,
                                                              const std::array<double,3> invariants);
  
  //! compute the elasticity tensor C = 2*dS/dC. Due to hyperelasticity there are symmetries C_{ijrs} = C_{jirs} and C_{ijrs} = C_{rsij} that leave 21 independent values.
  ElasticityTensor computeElasticityTensor(const double pressure, 
                                           const double pressureTilde,
                                           const std::array<Vec3,3> &rightCauchyGreen,
                                           const std::array<Vec3,3> &inverseRightCauchyGreen,
                                           const std::array<Vec3,3> &fictitiousPK2Stress,
                                           const std::array<Vec3,3> &pk2StressIsochoric,
                                           const double deformationGradientDeterminant,
                                           const std::array<double,2> reducedInvariants);
  
  //! compute the pressure from displacements, using the formula 2.8 of Sussman and Bathe "A finite element formulation for nonlinear incompressible elastic and inelastic analysis"
  double computePressureFromDisplacements(double deformationGradientDeterminant, 
                                          const std::array<Vec3,3> &rightCauchyGreen, 
                                          const std::array<Vec3,3> &PK2Stress);
};
 
};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.tpp"
