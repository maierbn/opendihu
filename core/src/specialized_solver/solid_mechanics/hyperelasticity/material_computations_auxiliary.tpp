#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include <Vc/Vc>

#include "equation/mooney_rivlin_incompressible.h"
#include "utility/math_utility.h"

namespace SpatialDiscretization
{

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeDeformationGradient(const std::array<VecD<3,double_v_t>,DisplacementsFunctionSpace::nDofsPerElement()> &displacements,
                           const Tensor2<3,double_v_t> &inverseJacobianMaterial,
                           const std::array<double, 3> xi
                          )
{
  // compute the deformation gradient x_i,j = δ_ij + u_i,j
  // where j is dimensionColumn and i is component of the used Vec3's


  VLOG(3) << "compute deformation gradient Fij, displacements: " << displacements;

  const int D = 3;
  Tensor2<D,double_v_t> deformationGradient;

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // compute du_i/dX_dimensionColumn for all i at once
    VecD<3,double_v_t> du_dX({0});    // vector contains entries for all i, du_dXj with j = dimensionColumn

    VLOG(3) << " j = " << dimensionColumn;

    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      VLOG(3) << "  M = " << dofIndex;

      // compute dphi_dofIndex/dX_dimensionColumn
      double_v_t dphi_dX = 0;
      for (int l = 0; l < 3; l++)
      {
        VLOG(3) << "   l = " << l;
        double_v_t dphi_dxil = DisplacementsFunctionSpace::dphi_dxi(dofIndex, l, xi);
        double_v_t dxil_dX = inverseJacobianMaterial[dimensionColumn][l];     // inverseJacobianMaterial[j][l] = J_lj = dxi_l/dX_j

        VLOG(3) << "     dphi_dxil = " << dphi_dxil << ", dxil_dX = " << dxil_dX;

        // multiply dphi/dxi with dxi/dX to obtain dphi/dX
        dphi_dX += dphi_dxil * dxil_dX;
      }
      VLOG(3) << "   dphi_dX = " << dphi_dX;

      VLOG(3) << "   displ_M: " << displacements[dofIndex];
      du_dX += dphi_dX * displacements[dofIndex];   // vector-valued addition
    }
    VLOG(3) << " du_dXj: " << du_dX;

    deformationGradient[dimensionColumn] = du_dX;

    // add Kronecker delta to obtain x_i,j = delta_ij + u_i,j
    deformationGradient[dimensionColumn][dimensionColumn] += 1;
  }

  VLOG(3) << "deformationGradient: " << deformationGradient;
  return deformationGradient;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeDeformationGradientTimeDerivative(const std::array<VecD<3,double_v_t>,DisplacementsFunctionSpace::nDofsPerElement()> &velocities,
                                         const Tensor2<3,double_v_t> &inverseJacobianMaterial,
                                         const std::array<double, 3> xi
                                        )
{
  // compute the time derivative of the deformation gradient d(x_i,j)/dt = d/dt(δ_ij + u_i,j) = v_i,j
  // where j is dimensionColumn and i is component of the used Vec3's


  VLOG(3) << "compute deformation gradient time derivative Fdot_ij, velocities: " << velocities;

  const int D = 3;
  Tensor2<D,double_v_t> deformationGradientTimeDerivative;

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // compute du_i/dX_dimensionColumn for all i at once
    VecD<3,double_v_t> du_dX({0});    // vector contains entries for all i, du_dXj with j = dimensionColumn

    VLOG(3) << " j = " << dimensionColumn;

    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      VLOG(3) << "  M = " << dofIndex;

      // compute dphi_dofIndex/dX_dimensionColumn
      double_v_t dphi_dX = 0;
      for (int l = 0; l < 3; l++)
      {
        VLOG(3) << "   l = " << l;
        double_v_t dphi_dxil = DisplacementsFunctionSpace::dphi_dxi(dofIndex, l, xi);
        double_v_t dxil_dX = inverseJacobianMaterial[dimensionColumn][l];     // inverseJacobianMaterial[j][l] = J_lj = dxi_l/dX_j

        VLOG(3) << "     dphi_dxil = " << dphi_dxil << ", dxil_dX = " << dxil_dX;

        // multiply dphi/dxi with dxi/dX to obtain dphi/dX
        dphi_dX += dphi_dxil * dxil_dX;
      }
      VLOG(3) << "   dphi_dX = " << dphi_dX;

      VLOG(3) << "   displ_M: " << velocities[dofIndex];
      du_dX += dphi_dX * velocities[dofIndex];   // vector-valued addition
    }
    VLOG(3) << " du_dXj: " << du_dX;

    deformationGradientTimeDerivative[dimensionColumn] = du_dX;
  }

  VLOG(3) << "deformationGradientTimeDerivative: " << deformationGradientTimeDerivative;
  return deformationGradientTimeDerivative;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeRightCauchyGreenTensor(const Tensor2<3,double_v_t> &deformationGradient)
{
  // compute C = F^T*F where F is the deformationGradient and C is the right Cauchy-Green Tensor
  // the quantities are 3x3 tensors for the 3D case
  Tensor2<3,double_v_t> rightCauchyGreenTensor({std::array<double_v_t,3>({0})});

  // C_ji
  // loop over dimension, i.e. columns of right cauchy green tensor, i
  for (int dimensionColumn = 0; dimensionColumn < 3; dimensionColumn++)
  {
    // loop over row of tensor, j
    for (int dimensionRow = 0; dimensionRow < 3; dimensionRow++)
    {
      for (int k = 0; k < 3; k++)
      {
        // C_cr = C_rc += (F^T)_ck * F_kr = F_kc * F_kr
        rightCauchyGreenTensor[dimensionColumn][dimensionRow] +=
          deformationGradient[dimensionColumn][k] * deformationGradient[dimensionRow][k];
      }
    }
  }

  return rightCauchyGreenTensor;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
std::array<double_v_t,5> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeInvariants(const Tensor2<3,double_v_t> &rightCauchyGreen, const double_v_t rightCauchyGreenDeterminant, const VecD<3,double_v_t> fiberDirection)
{
  std::array<double_v_t,5> invariants;

  // I1 = tr(C)
  invariants[0] = 0.0;
  for (int i = 0; i < 3; i++)
  {
    invariants[0] += rightCauchyGreen[i][i];
  }

  // tr(C^2)
  double_v_t traceC2 = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      traceC2 += rightCauchyGreen[i][j] * rightCauchyGreen[j][i];
    }
  }

  // I2 = 1/2 * (tr(C)^2 - tr(C^2))
  invariants[1] = 0.5 * (MathUtility::sqr(invariants[0]) - traceC2);

  // I3 = det(C)
  invariants[2] = rightCauchyGreenDeterminant;

  // for a transversely isotropic material that also has 4th and 5th invariant
  if (Term::usesFiberDirection)
  {
    // I4 = a0 • C a0;
    double_v_t a0Ca0 = 0;
    for (int i = 0; i < 3; i++)
    {
      double_v_t ca0_i = 0;
      for (int j = 0; j < 3; j++)
      {
        ca0_i += rightCauchyGreen[j][i] * fiberDirection[j];
      }
      a0Ca0 += fiberDirection[i] * ca0_i;
    }

    invariants[3] = a0Ca0;

    // I5 = a0 • C^2 a0;

    double_v_t a0C2a0 = 0;
    for (int i = 0; i < 3; i++)
    {
      double_v_t c2a0_i = 0;
      for (int j = 0; j < 3; j++)
      {
        // compute C^2
        double_v_t c2_ij = 0;
        for (int k = 0; k < 3; k++)
        {
          // C^2_ij = C_ik * C_kj
          c2_ij += rightCauchyGreen[k][i] * rightCauchyGreen[j][k];
        }
        c2a0_i += c2_ij * fiberDirection[j];
      }
      a0C2a0 += fiberDirection[i] * c2a0_i;
    }

    invariants[4] = a0C2a0;
    //LOG(DEBUG) << "computed I4: " << invariants[3] << ", I5: " << invariants[4] << ", a0: " << fiberDirection;
  }

  return invariants;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
std::array<double_v_t,5> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeReducedInvariants(const std::array<double_v_t,5> invariants, const double_v_t deformationGradientDeterminant)
{
  std::array<double_v_t,5> reducedInvariants;

  // 3D: Fbar = J^{-1/3}*F such that det(Fbar) = 1
  // 2D: Fbar = J^{-1/2}*F such that det(Fbar) = 1

  // 3D: Cbar = J^{-2/3}*C such that det(Cbar) = 1
  // 2D: Cbar = J^{-2/2}*C = J^{-1}*C = J^{-1}*F^T*F = J^{-1}*J^{1/2}*J^{1/2}*Fbar^T*Fbar = Fbar^T*Fbar such that det(Cbar) = 1

  double factor23 = -2./3;
  double factor43 = -4./3;

  // for 2D, not used here
#if 0
  factor23 = -1./2;
  factor43 = -2./2;
#endif

  reducedInvariants[0] = MathUtility::pow(deformationGradientDeterminant, factor23) * invariants[0];
  reducedInvariants[1] = MathUtility::pow(deformationGradientDeterminant, factor43) * invariants[1];

  // for a transversely isotropic material that also has 4th and 5th invariant
  if (Term::usesFiberDirection)
  {
    reducedInvariants[2] = 1;  // not used, because I3 = det C = 1 constant
    reducedInvariants[3] = MathUtility::pow(deformationGradientDeterminant, factor23) * invariants[3];
    reducedInvariants[4] = MathUtility::pow(deformationGradientDeterminant, factor43) * invariants[4];
  }

  if (Vc::any_of(deformationGradientDeterminant <= 0))
  {
    LOG(ERROR) << "J=det F is negative: " << deformationGradientDeterminant << ". Result will be unphysical.\n"
      << "For dynamic problems, reduce time step width, for static problems, add smaller \"loadFactors\" or reduce load.";

    Vc::where(deformationGradientDeterminant <= 0) | reducedInvariants[0] = 3;
    Vc::where(deformationGradientDeterminant <= 0) | reducedInvariants[1] = 0;
    Vc::where(deformationGradientDeterminant <= 0) | reducedInvariants[3] = 3;
    Vc::where(deformationGradientDeterminant <= 0) | reducedInvariants[4] = 0;
  }

  return reducedInvariants;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computePK2Stress(double_v_t &pressure,                                   //< [in] pressure value p
                 const Tensor2<3,double_v_t> &rightCauchyGreen,         //< [in] C
                 const Tensor2<3,double_v_t> &inverseRightCauchyGreen,  //< [in] C^{-1}
                 const std::array<double_v_t,5> invariants,             //< [in] the strain invariants I_1, ..., I_5
                 const std::array<double_v_t,5> reducedInvariants,      //< [in] the reduced invariants Ibar_1, ..., Ibar_5
                 const double_v_t deformationGradientDeterminant,       //< [in] J = det(F)
                 VecD<3,double_v_t> fiberDirection,                     //< [in] a0, direction of fibers
                 Tensor2<3,double_v_t> &fictitiousPK2Stress,            //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                 Tensor2<3,double_v_t> &pk2StressIsochoric              //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                )
{
  // compute the PK2 stress tensor as S=2*dPsi/dC
  // for explanation see pdf document

  // if the alternative coupled form of the strain energy function, ψ(C), strainEnergyDensityFunctionCoupledDependentOnC, is considered
  const bool usesFormulationWithC = typeid(decltype(Term::strainEnergyDensityFunctionCoupledDependentOnC)) != typeid(decltype(INT(0)));

  // Ibar1, Ibar2, Ibar4, Ibar5, J, I1, I2, I3, C11, C12, C13, C22, C23, C33
  std::vector<double_v_t> parameterVector = {
    reducedInvariants[0], reducedInvariants[1], reducedInvariants[3], reducedInvariants[4],  // Ibar1, Ibar2, Ibar4, Ibar5
    deformationGradientDeterminant,                                                          // J
    invariants[0], invariants[1], invariants[2],                                             // I1, I2, I3
    rightCauchyGreen[0][0], rightCauchyGreen[1][0], rightCauchyGreen[2][0],                  // C11, C12, C13
    rightCauchyGreen[1][1], rightCauchyGreen[2][1], rightCauchyGreen[2][2]                   // C22, C23, C33
  };

  // reduced invariants, arguments of `strainEnergyDensityFunctionIsochoric`
  // compute factors for decoupled form
  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);

  const double_v_t Ibar1 = reducedInvariants[0];

  const double_v_t dPsi_dIbar1 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar1Expression, parameterVector);
  const double_v_t dPsi_dIbar2 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar2Expression, parameterVector);

  const double_v_t J = deformationGradientDeterminant;

  double_v_t decoupledFormFactor1 = 2*(dPsi_dIbar1 + Ibar1*dPsi_dIbar2);
  double_v_t decoupledFormFactor2 = -2*dPsi_dIbar2;

  // terms for 4th and 5th invariants
  double_v_t decoupledFormFactor4 = 0;
  double_v_t decoupledFormFactor5 = 0;

  if (Term::usesFiberDirection)
  {
    auto dPsi_dIbar4Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar4);
    auto dPsi_dIbar5Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar5);

    const double_v_t dPsi_dIbar4 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar4Expression, parameterVector);
    const double_v_t dPsi_dIbar5 = ExpressionHelper<double_v_t>::apply(dPsi_dIbar5Expression, parameterVector);

    decoupledFormFactor4 = 2*dPsi_dIbar4;
    decoupledFormFactor5 = 2*dPsi_dIbar5;
  }

  // compute factors for coupled form
  const double_v_t I1 = invariants[0];
  const double_v_t I3 = invariants[2];

  auto dPsi_dI1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupled, Term::I1);
  auto dPsi_dI2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupled, Term::I2);
  auto dPsi_dI3Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupled, Term::I3);
  const double_v_t dPsi_dI1 = ExpressionHelper<double_v_t>::apply(dPsi_dI1Expression, parameterVector);
  const double_v_t dPsi_dI2 = ExpressionHelper<double_v_t>::apply(dPsi_dI2Expression, parameterVector);
  const double_v_t dPsi_dI3 = ExpressionHelper<double_v_t>::apply(dPsi_dI3Expression, parameterVector);

  VLOG(2) << "coupled term: " << Term::strainEnergyDensityFunctionCoupled;
  VLOG(2) << "invariants: I1: " << Term::I1 << " = " << I1 << ", I2: " << Term::I2 << ", I3: " << Term::I3 << " = " << I3;
  VLOG(2) << "∂ψ/∂I1: " << dPsi_dI1Expression << " = " << dPsi_dI1;
  VLOG(2) << "∂ψ/∂I2: " << dPsi_dI2Expression << " = " << dPsi_dI2;
  VLOG(2) << "∂ψ/∂I3: " << dPsi_dI3Expression << " = " << dPsi_dI3;

  double_v_t coupledFormFactor1 = 2*(dPsi_dI1 + I1 * dPsi_dI2);
  double_v_t coupledFormFactor2 = -2*dPsi_dI2;
  double_v_t coupledFormFactor3 = 2 * I3 * dPsi_dI3;


  // for compressible material, the stress, p, is given by a constitutive equation
  // for incompressible material, the stress, p, is an unknown that will be solved for
  if (!Term::isIncompressible)
  {
    auto dPsi_dJExpression = SEMT::deriv_t(Term::strainEnergyDensityFunctionVolumetric, Term::J);
    const double_v_t dPsi_dJ = ExpressionHelper<double_v_t>::apply(dPsi_dJExpression, parameterVector);
    VLOG(1) << "J = " << Term::J << " = " << J << ", Psi_vol: " << Term::strainEnergyDensityFunctionVolumetric << ", dΨ/dJ: " << dPsi_dJExpression << " = " << dPsi_dJ;
    pressure = dPsi_dJ;
  }

  double factor23 = -2./3;
  double_v_t factorJ23 = MathUtility::pow(J, factor23);

  Tensor2<3,double_v_t> pK2Stress;

  // compute fictitiousPK2Stress:
  // Sbar = factor1*I + factor2*Cbar, Cbar = J^{-2/3}*C

  // row index
  for (int i=0; i<3; i++)
  {
    // column index
    for (int j=0; j<3; j++)
    {
      int delta_ij = (i == j? 1 : 0);
      const double_v_t cBar = factorJ23 * rightCauchyGreen[j][i];
      fictitiousPK2Stress[j][i] = decoupledFormFactor1 * delta_ij + decoupledFormFactor2 * cBar;

      //if (i == j)
      //LOG(DEBUG) << "  C_" << i << j << " = " << rightCauchyGreen[j][i] << ", J=" << J << ", factorJ23 = " << factorJ23 << ", SBar_" << i << j << ": " << fictitiousPK2Stress[j][i]
      //     << " = " << factor1 << " * " << delta_ij << " + " << factor2 << "*" << cBar;
    }
  }

  //LOG(DEBUG) << "in computePK2Stress, factor1: " << factor1 << ", factor2: " << factor2 << ", C: " << rightCauchyGreen << ", Sbar: " << fictitiousPK2Stress << ", reducedInvariants: " << reducedInvariants;

  // add terms for 4th and 5th invariants
  if (Term::usesFiberDirection)
  {
    // add term for 4th invariant
    for (int i=0; i<3; i++)     // alternative indices: A
    {
      // column index
      for (int j=0; j<3; j++)     // alternative indices: B
      {
        fictitiousPK2Stress[j][i] += decoupledFormFactor4 * fiberDirection[i] * fiberDirection[j];
      }
    }

    //LOG(DEBUG) << "fiberDirection: " << fiberDirection << ", decoupledFormFactor4: " << decoupledFormFactor4
    //  << ", dPsi_dIbar4Expression: " << dPsi_dIbar4Expression << ", invariants: " << reducedInvariantsVector;

    // add term for 5th invariant
    for (int i=0; i<3; i++)     // alternative indices: A
    {
      // compute a0 Cbar, a0cbar_i = a0_k*C_ki
      double_v_t cbara0_i = 0;
      for (int k=0; k<3; k++)
      {
        cbara0_i += fiberDirection[k] * rightCauchyGreen[i][k];
      }

      // column index
      for (int j=0; j<3; j++)     // alternative indices: B
      {
        // compute Cbar a0, cbara0_j = C_jk * a0_k
        double_v_t cbara0_j = 0;
        for (int k=0; k<3; k++)
        {
          cbara0_j += rightCauchyGreen[k][j] * fiberDirection[k];
        }

        // Sbar += gamma_5 * (a0 dyad Cbar a0  +  a0 Cbar dyad a0)
        fictitiousPK2Stress[j][i] += decoupledFormFactor5 * (fiberDirection[i] * cbara0_j + cbara0_i * fiberDirection[j]);
      }
    }
  }

  // decoupled form of strain energy function
  // Holzapfel p.234
  // S = S_vol + S_iso
  // S_vol = J*p*C^-1
  // S_iso = 2*dPsi_iso/dC = J^(-2/3) P : Sbar   (P = II - 1/3(C^{-1} dyad C): Holzapfel p.229)
  // pSbar = P : Sbar, P: 4th order tensor, Sbar: 2nd order tensor, note: A:B = A_ijkl*B_kl*e_i dyad e_j

  // compute S = S_vol + S_iso
  // row index
  for (int i=0; i<3; i++)     // alternative indices: A
  {
    // column index
    for (int j=0; j<3; j++)     // alternative indices: B
    {
      // volumetric stress
      const double_v_t sVol = J * pressure * inverseRightCauchyGreen[j][i];         // S_vol = J * p * C^{-1}_AB

      // compute P : Sbar
      double_v_t pSbar = 0;
      // row index
      for (int k=0; k<3; k++)        // alternative indices: C
      {
        const int delta_ik = (i == k? 1 : 0);

        // column index
        for (int l=0; l<3; l++)            // alternative indices: D
        {
          const int delta_jl = (j == l? 1 : 0);

          // this is a non-symmetric version for Ii but it is also correct, the symmetric version would be given by Ii = (δ_AC*δ_BD + δ_AD*δ_BC) / 2
          const int Ii = delta_ik * delta_jl;
          const double_v_t Cc = inverseRightCauchyGreen[j][i] * rightCauchyGreen[l][k];     // CC = C^{-1}_AB * C_CD
          const double_v_t Pp = (Ii - 1./3 * Cc);

          //LOG(DEBUG) << "    PP_" << i << j << k << l << " = " << Pp << " = " << Ii << "-1/3*" << Cc << " (" << inverseRightCauchyGreen[j][i] << "," << rightCauchyGreen[l][k] << "), Ii: " << Ii;

          pSbar += Pp * fictitiousPK2Stress[l][k];
        }
      }


      // isochoric stress
      const double_v_t sIso = factorJ23 * pSbar;
      pk2StressIsochoric[j][i] = sIso;

      //VLOG(2) << "    PSbar_" << i << j << ": " << pSbar << ", J^-2/3: " << factorJ23 << ", Siso_" << i << j << ": " << sIso;
      VLOG(2) << "   sVol: " << sVol << ", sIso: " << sIso;

      // total stress is sum of volumetric and isochoric part, sVol = J*p*C^{-1}_AB, sIso = j^{-2/3}*(Ii-1/3*Cc)*Sbar
      pK2Stress[j][i] = sVol + sIso;

      // debugging
      //pK2Stress[j][i] = sIso;

      //VLOG(2) << "set pk2Stress_" << i << j << " = " << pK2Stress[j][i];

      //if (i == j)
      //  LOG(DEBUG) << "  ccs: " << ccs << " C:Sbar: " << cSbar << ", factorJ23: " << factorJ23 << ", Svol_" << i << j << " = " << sVol << ", Siso_" << i << j << " = " << sIso << ", S = " << pK2Stress[j][i];
    }  // j
  }  // i

  // coupled form of strain energy function
  // compute S += γ1*I + γ2*C + γ3*C^-1   (Holzapfel, p.248)
  // row index
  for (int i=0; i<3; i++)     // alternative indices: A
  {
    // column index
    for (int j=0; j<3; j++)     // alternative indices: B
    {
      const int delta_ij = (i == j? 1 : 0);

      pK2Stress[j][i] += coupledFormFactor1 * delta_ij + coupledFormFactor2 * rightCauchyGreen[j][i] + coupledFormFactor3 * inverseRightCauchyGreen[j][i];
    }
  }
  VLOG(2) << "   δ1: " << coupledFormFactor1 << ", δ2: " << coupledFormFactor2 << ", δ3: " << coupledFormFactor3 << ", after adding coupled form terms: " << pK2Stress;

  // if the formulation includes a Ψ(C) term, we need to add S = 2 * ∂Ψ(C)/∂C
  if (usesFormulationWithC)
  {
    // C = F^T F is symmetric, we only use the entries C11, C12, C13, C22, C23, C33
    auto dPsi_dC11Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C11);
    auto dPsi_dC12Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C12);
    auto dPsi_dC13Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C13);
    auto dPsi_dC22Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C22);
    auto dPsi_dC23Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C23);
    auto dPsi_dC33Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C33);

    const double_v_t dPsi_dC11 = ExpressionHelper<double_v_t>::apply(dPsi_dC11Expression, parameterVector);
    const double_v_t dPsi_dC12 = ExpressionHelper<double_v_t>::apply(dPsi_dC12Expression, parameterVector);
    const double_v_t dPsi_dC13 = ExpressionHelper<double_v_t>::apply(dPsi_dC13Expression, parameterVector);
    const double_v_t dPsi_dC22 = ExpressionHelper<double_v_t>::apply(dPsi_dC22Expression, parameterVector);
    const double_v_t dPsi_dC23 = ExpressionHelper<double_v_t>::apply(dPsi_dC23Expression, parameterVector);
    const double_v_t dPsi_dC33 = ExpressionHelper<double_v_t>::apply(dPsi_dC33Expression, parameterVector);

    VLOG(2) << "formulation with C, add S=(" << dPsi_dC11 << "," << 2*dPsi_dC12 << "," << 2*dPsi_dC13 << "," << 2*dPsi_dC22 << "," << 2*dPsi_dC23 << "," << 2*dPsi_dC33 << ")";

    // add contribution to S = 2 * ∂Ψ(C)/∂C
    pK2Stress[0][0] += 2*dPsi_dC11;   // S11
    pK2Stress[1][0] += 2*dPsi_dC12;   // S12
    pK2Stress[2][0] += 2*dPsi_dC13;   // S13
    pK2Stress[0][1] += 2*dPsi_dC12;   // S21 = S12
    pK2Stress[1][1] += 2*dPsi_dC22;   // S22
    pK2Stress[2][1] += 2*dPsi_dC23;   // S23
    pK2Stress[2][0] += 2*dPsi_dC13;   // S31 = S13
    pK2Stress[2][1] += 2*dPsi_dC23;   // S32 = S23
    pK2Stress[2][2] += 2*dPsi_dC33;   // S33
  }

  //for debugging, check symmetry of PK2 stress and if it is correct according to Mooney-Rivlin formula
  if (VLOG_IS_ON(2))
  {
    const double_v_t errorTolerance = 1e-14;
    bool pK2IsSymmetric = true;
    for (int a=0; a<3; a++)
    {
      for (int b=0; b<3; b++)
      {
        if (Vc::any_of(MathUtility::abs(pK2Stress[a][b] - pK2Stress[b][a]) > errorTolerance))
        {
          LOG(ERROR) << "pK2Stress[" <<a<< "][" <<b<< "] != pK2Stress[" <<b<< "][" <<a<< "] (" <<pK2Stress[b][a]<< " != " <<pK2Stress[a][b]<< ")";
          pK2IsSymmetric = false;
        }
      }
    }
    if (pK2IsSymmetric)
      VLOG(2) << "PK2 stress tensor is symmetric!";

    // check if PK2 is the same as from explicit formula for Mooney-Rivlin (p.249)

    // explicit formula in Holzapfel p.249

    double factor23 = -2./3;
    const double_v_t factorJ23 = MathUtility::pow(deformationGradientDeterminant, factor23);

    const double c0 = SEMT::Parameter<0>::get_value();
    const double c1 = SEMT::Parameter<1>::get_value();
    //const double c0 = PARAM(0).get_value();   //< material parameter
    //const double c1 = PARAM(1).get_value();   //< material parameter

    const double_v_t Ibar1 = reducedInvariants[0];

    const double_v_t gamma1 = 2*(c0 + c1*Ibar1);
    const double_v_t gamma2 = -2*c1;

    bool mismatch = false;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        const int delta_ab = (a == b? 1 : 0);
        double_v_t cBar = factorJ23 * rightCauchyGreen[b][a];
        double_v_t sBar = gamma1*delta_ab + gamma2*cBar;

        if (Vc::any_of(MathUtility::abs(sBar - fictitiousPK2Stress[b][a]) > errorTolerance))
        {
          LOG(ERROR) << "mismatch in Sbar_" << a << b << ": derived: " << fictitiousPK2Stress << ", Mooney Rivlin explicit formula: " << sBar;
          mismatch = true;
        }
      }
    }
    if (!mismatch)
      LOG_N_TIMES(2,DEBUG) << "Sbar is correct!";
  }

#ifndef NDEBUG
  if (MathUtility::containsNanOrInf(pK2Stress))
  {
    LOG(FATAL) << "PK2stress contains nan: " << pK2Stress;
  }
#endif

  return pK2Stress;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computePK2StressField()
{
  //LOG(TRACE) << "computePK2StressField";

  //this->data_.pK2Stress()->startGhostManipulation();
  this->data_.pK2Stress()->zeroGhostBuffer();
  this->data_.materialTraction()->zeroGhostBuffer();
  this->data_.materialTraction()->zeroEntries();

  this->data_.deformationGradient()->zeroGhostBuffer();
  this->data_.deformationGradientTimeDerivative()->zeroGhostBuffer();

  this->data_.geometryReference()->setRepresentationGlobal();
  this->data_.geometryReference()->startGhostManipulation();

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nPressureDofsPerElement = PressureFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();

  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    // get geometry field of reference configuration
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference()->getElementValues(elementNoLocalv, geometryReferenceValues);

    // get displacements field values for element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> displacementsValues;
    this->data_.displacements()->getElementValues(elementNoLocalv, displacementsValues);

    // get displacements field values for element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> velocitiesValues;
    this->data_.velocities()->getElementValues(elementNoLocalv, velocitiesValues);

    std::array<double_v_t,nPressureDofsPerElement> pressureValuesCurrentElement;
    this->data_.pressure()->getElementValues(elementNoLocalv, pressureValuesCurrentElement);

    std::array<Vec3_v_t,nDisplacementsDofsPerElement> elementalDirectionValues;
    this->data_.fiberDirection()->getElementValues(elementNoLocalv, elementalDirectionValues);

    // get indices of element-local dofs
    std::array<dof_no_v_t,27> dofNosLocal = this->displacementsFunctionSpace_->getElementDofNosLocal(elementNoLocalv);

    //LOG(DEBUG) << "el " << elementNoLocal << ", geometryRef: " << geometryReferenceValues << ", displacements: " << displacementsValues << ", p: " << pressureValuesCurrentElement;

    // loop over nodes of this element
    for (int elementalNodeNo = 0; elementalNodeNo < 27; elementalNodeNo++)
    {
      dof_no_v_t dofNoLocal = dofNosLocal[elementalNodeNo];

      // get parameter values of current sampling point
      int indexX = elementalNodeNo % 3;
      int indexY = int((elementalNodeNo % 9) / 3);
      int indexZ = int(elementalNodeNo / 9);
      Vec3 xi({});

      switch (indexX)
      {
      case 0:
        xi[0] = 0.0;
        break;
      case 1:
        xi[0] = 0.5;
        break;
      case 2:
        xi[0] = 1.0;
        break;
      };
      switch (indexY)
      {
      case 0:
        xi[1] = 0.0;
        break;
      case 1:
        xi[1] = 0.5;
        break;
      case 2:
        xi[1] = 1.0;
        break;
      };
      switch (indexZ)
      {
      case 0:
        xi[2] = 0.0;
        break;
      case 1:
        xi[2] = 0.5;
        break;
      case 2:
        xi[2] = 1.0;
        break;
      };

      Tensor2_v_t<D> jacobianMaterial;
      double_v_t jacobianDeterminant;
      Tensor2_v_t<D> inverseJacobianMaterial;
      Tensor2_v_t<D> deformationGradient;
      double_v_t deformationGradientDeterminant;

      for (int nTries = 0; nTries < 3; nTries++)
      {
        // compute the 3x3 jacobian of the parameter space to world space mapping
        jacobianMaterial = DisplacementsFunctionSpace::computeJacobian(geometryReferenceValues, xi);
        inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);

        // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
        // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

        // F
        deformationGradient = this->computeDeformationGradient(displacementsValues, inverseJacobianMaterial, xi);
        deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J

        if (Vc::all_of(deformationGradientDeterminant > 0.2))
          break;

        // if J=det(F) is negative, move the point in the element, xi, a bit more to the center
        for (int i = 0; i < 3; i++)
        {
          xi[i] = 0.5 + (xi[0]-0.5)*0.9;
        }

#ifndef NDEBUG
        LOG(DEBUG) << "element " << elementNoLocal << ", J=" << deformationGradientDeterminant << "," << jacobianDeterminant
          << ", displacementsValues: " << displacementsValues[0] << "," << displacementsValues[1]
          << ", deformationGradient: " << deformationGradient << ", inverseJacobianMaterial: " << inverseJacobianMaterial
          << ", geometryReferenceValues: " << geometryReferenceValues[0] << "," << geometryReferenceValues[1]
          << ", " << nTries << " retry with xi=" << xi;
#endif
      }

      if (Vc::any_of(deformationGradientDeterminant < 0))
      {
        LOG(ERROR) << "J = det(F) = " << deformationGradientDeterminant << " is negative, in computation of PK2 stresses.\n"
          << "Element no. " << elementNoLocal << ", xi=" << xi << ", det(material jacobian): " << jacobianDeterminant
          << ", displacementsValues: " << displacementsValues[0] << "," << displacementsValues[1]
          << ", deformationGradient: " << deformationGradient << ", inverseJacobianMaterial: " << inverseJacobianMaterial
          << ", geometryReferenceValues: " << geometryReferenceValues[0] << "," << geometryReferenceValues[1];
      }

      // compute Fdot values
      Tensor2_v_t<D> Fdot = computeDeformationGradientTimeDerivative(velocitiesValues, inverseJacobianMaterial, xi);

      // store F values
      std::array<double_v_t,9> deformationGradientValues;
      std::array<double_v_t,9> deformationGradientTimeDerivativeValues;
      for (int j = 0; j < 3; j++)
      {
        for (int i = 0; i < 3; i++)
        {
          // row-major
          deformationGradientValues[j*3+i] = deformationGradient[j][i];
          deformationGradientTimeDerivativeValues[j*3+i] = Fdot[j][i];
        }
      }
      this->data_.deformationGradient()->setValue(dofNoLocal, deformationGradientValues, INSERT_VALUES);
      this->data_.deformationGradientTimeDerivative()->setValue(dofNoLocal, deformationGradientTimeDerivativeValues, INSERT_VALUES);


      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

      // invariants
      std::array<double_v_t,5> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant, fiberDirection);  // I_1, I_2, I_3, I_4, I_5
      std::array<double_v_t,5> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2, Ibar_4, Ibar_5

      // pressure is the separately interpolated pressure for mixed formulation
      double_v_t pressure = pressureFunctionSpace->interpolateValueInElement(pressureValuesCurrentElement, xi);

      // checking for nans in debug mode
#ifndef NDEBUG
      if (MathUtility::containsNanOrInf(inverseJacobianMaterial))
        LOG(FATAL) << "inverseJacobianMaterial contains nan: " << inverseJacobianMaterial << ", jacobianMaterial: " << jacobianMaterial;

      if (MathUtility::containsNanOrInf(deformationGradient))
        LOG(FATAL) << "deformationGradient contains nan: " << deformationGradient;

      if (MathUtility::containsNanOrInf(rightCauchyGreen))
        LOG(FATAL) << "rightCauchyGreen contains nan: " << rightCauchyGreen;

      if (MathUtility::containsNanOrInf(inverseRightCauchyGreen))
        LOG(FATAL) << "inverseRightCauchyGreen contains nan: " << inverseRightCauchyGreen;

      if (MathUtility::containsNanOrInf(invariants))
        LOG(FATAL) << "invariants contains nan: " << invariants;

      if (MathUtility::containsNanOrInf(deformationGradientDeterminant))
        LOG(FATAL) << "deformationGradientDeterminant contains nan: " << deformationGradientDeterminant;

      if (MathUtility::containsNanOrInf(reducedInvariants))
        LOG(FATAL) << "reducedInvariants contains nan: " << reducedInvariants << ", deformationGradient: " << deformationGradient << ", deformationGradientDeterminant: " << deformationGradientDeterminant;

      if (MathUtility::containsNanOrInf(pressure))
        LOG(FATAL) << "pressure contains nan: " << pressure;
#endif

      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2_v_t<D> fictitiousPK2Stress;   // Sbar
      Tensor2_v_t<D> pk2StressIsochoric;    // S_iso
      Tensor2_v_t<D> pK2Stress = this->computePK2Stress(pressure, rightCauchyGreen, inverseRightCauchyGreen, invariants, reducedInvariants, deformationGradientDeterminant, fiberDirection,
                                                      fictitiousPK2Stress, pk2StressIsochoric
                                                    );

      std::array<double_v_t,6> valuesInVoigtNotation({pK2Stress[0][0], pK2Stress[1][1], pK2Stress[2][2], pK2Stress[0][1], pK2Stress[1][2], pK2Stress[0][2]});

      //LOG(INFO) << "node " << dofNoLocal << " pk2: " << valuesInVoigtNotation;
      this->data_.pK2Stress()->setValue(dofNoLocal, valuesInVoigtNotation, INSERT_VALUES);

      // compute surface traction
      // get normal
      if (indexZ == 0)
      {
        // bottom node
        Vec3_v_t normal = displacementsFunctionSpace->getNormal(Mesh::face_t::face2Minus, elementNoLocalv, xi);

        // compute traction by Cauchy theorem T = S n
        Vec3_v_t traction = pK2Stress * normal;

        // set value in material traction
        this->data_.materialTraction()->setValue(dofNoLocal, traction, INSERT_VALUES);

      }
      else if (indexZ == 2)
      {
        // top node
        Vec3_v_t normal = displacementsFunctionSpace->getNormal(Mesh::face_t::face2Plus, elementNoLocalv, xi);

        // compute traction by Cauchy theorem T = S n
        Vec3_v_t traction = pK2Stress * normal;

        // set value in material traction
        this->data_.materialTraction()->setValue(dofNoLocal, traction, INSERT_VALUES);
      }
    }
  }
  this->data_.pK2Stress()->zeroGhostBuffer();
  this->data_.pK2Stress()->finishGhostManipulation();
  this->data_.pK2Stress()->startGhostManipulation();

  this->data_.materialTraction()->zeroGhostBuffer();
  this->data_.materialTraction()->finishGhostManipulation();
  this->data_.materialTraction()->startGhostManipulation();

  this->data_.deformationGradient()->zeroGhostBuffer();
  this->data_.deformationGradient()->finishGhostManipulation();
  this->data_.deformationGradient()->startGhostManipulation();

  this->data_.deformationGradientTimeDerivative()->zeroGhostBuffer();
  this->data_.deformationGradientTimeDerivative()->finishGhostManipulation();
  this->data_.deformationGradientTimeDerivative()->startGhostManipulation();
}

//! compute the material elasticity tensor
template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeElasticityTensor(const Tensor2<3,double_v_t> &rightCauchyGreen,         //< [in] C
                        const Tensor2<3,double_v_t> &inverseRightCauchyGreen,  //< [in] C^{-1}
                        double_v_t deformationGradientDeterminant,             //< [in] J = det(F)
                        double_v_t pressure,                                   //< [in] pressure value p
                        std::array<double_v_t,5> invariants,                   //< [in] the invariants I1, ..., I5
                        std::array<double_v_t,5> reducedInvariants,            //< [in] the reduced invariants Ibar_1, ..., Ibar_5
                        const Tensor2<3,double_v_t> &fictitiousPK2Stress,      //< [in] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                        const Tensor2<3,double_v_t> &pk2StressIsochoric,       //< [in] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                        VecD<3,double_v_t> fiberDirection,                     //< [in] a0, direction of fibers
                        Tensor4<3,double_v_t> &fictitiousElasticityTensor,     //< [out] fictitious Elasticity tensor CCbar_{ABCD}
                        Tensor4<3,double_v_t> &elasticityTensorIso,            //< [out] CCiso_{ABCD}
                        Tensor4<3,double_v_t> &elasticityTensor                //< [out] elasticity tensor CC_{ABCD}
                       )
{
  // compute the elasticity tensor as CC=2*dS(C)/dC
  // for explanation see pdf document
  const int D = 3;

  // if the alternative coupled form of the strain energy function, ψ(C), strainEnergyDensityFunctionCoupledDependentOnC, is considered
  const bool usesFormulationWithC = typeid(decltype(Term::strainEnergyDensityFunctionCoupledDependentOnC)) != typeid(decltype(INT(0)));

  // Ibar1, Ibar2, Ibar4, Ibar5, J, I1, I2, I3, C11, C12, C13, C22, C23, C33
  std::vector<double_v_t> parameterVector = {
    reducedInvariants[0], reducedInvariants[1], reducedInvariants[3], reducedInvariants[4],  // Ibar1, Ibar2, Ibar4, Ibar5
    deformationGradientDeterminant,                                                          // J
    invariants[0], invariants[1], invariants[2],                                             // I1, I2, I3
    rightCauchyGreen[0][0], rightCauchyGreen[1][0], rightCauchyGreen[2][0],                  // C11, C12, C13
    rightCauchyGreen[1][1], rightCauchyGreen[2][1], rightCauchyGreen[2][2]                   // C22, C23, C33
  };

  // compute preliminary variables that are independent of the indices a,b,c,d
  // decoupled form of strain energy function
  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);
  auto dPsi_dIbar5Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar5);

  auto d2Psi_dIbar1Ibar1Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar1);
  auto d2Psi_dIbar1Ibar2Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar2);
  auto d2Psi_dIbar2Ibar2Expression = SEMT::deriv_t(dPsi_dIbar2Expression, Term::Ibar2);

  const double_v_t Ibar1 = reducedInvariants[0];

  const double_v_t dPsi_dIbar2       = ExpressionHelper<double_v_t>::apply(dPsi_dIbar2Expression, parameterVector);
  const double_v_t dPsi_dIbar5       = ExpressionHelper<double_v_t>::apply(dPsi_dIbar5Expression, parameterVector);
  const double_v_t d2Psi_dIbar1Ibar1 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar1Expression, parameterVector);
  const double_v_t d2Psi_dIbar1Ibar2 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar2Expression, parameterVector);
  const double_v_t d2Psi_dIbar2Ibar2 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar2Ibar2Expression, parameterVector);

  const double_v_t J = deformationGradientDeterminant;

  const double_v_t factorJ23 = MathUtility::pow(J,-2./3);   // J^{-2/3}
  const double_v_t factorJ43 = MathUtility::pow(J,-4./3);   // J^{-4/3}

  // compute factors for Cbar
  const double_v_t decoupledFormFactor1 = 4*(d2Psi_dIbar1Ibar1 + 2*Ibar1*d2Psi_dIbar1Ibar2 + dPsi_dIbar2 + MathUtility::sqr(Ibar1)*d2Psi_dIbar2Ibar2);
  const double_v_t decoupledFormFactor2 = -4*(d2Psi_dIbar1Ibar2 + Ibar1*d2Psi_dIbar2Ibar2);
  const double_v_t decoupledFormFactor3 = 4*d2Psi_dIbar2Ibar2;
  const double_v_t decoupledFormFactor4 = -4*dPsi_dIbar2;

  // terms for 4th and 5th invariant
  double_v_t decoupledFormFactor5 = 0;
  double_v_t decoupledFormFactor6 = 0;
  double_v_t decoupledFormFactor7 = 0;
  double_v_t decoupledFormFactor8 = 0;
  double_v_t decoupledFormFactor9 = 0;
  double_v_t decoupledFormFactor10 = 0;
  double_v_t decoupledFormFactor11 = 0;
  double_v_t decoupledFormFactor12 = 0;

  if (Term::usesFiberDirection)
  {
    auto dPsi_dIbar4Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar4);

    auto d2Psi_dIbar1Ibar4Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar4);
    auto d2Psi_dIbar2Ibar4Expression = SEMT::deriv_t(dPsi_dIbar2Expression, Term::Ibar4);
    auto d2Psi_dIbar4Ibar4Expression = SEMT::deriv_t(dPsi_dIbar4Expression, Term::Ibar4);
    auto d2Psi_dIbar1Ibar5Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar5);
    auto d2Psi_dIbar2Ibar5Expression = SEMT::deriv_t(dPsi_dIbar2Expression, Term::Ibar5);
    auto d2Psi_dIbar5Ibar5Expression = SEMT::deriv_t(dPsi_dIbar5Expression, Term::Ibar5);
    auto d2Psi_dIbar4Ibar5Expression = SEMT::deriv_t(dPsi_dIbar4Expression, Term::Ibar5);

    const double_v_t d2Psi_dIbar1Ibar4 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar4Expression, parameterVector);
    const double_v_t d2Psi_dIbar2Ibar4 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar2Ibar4Expression, parameterVector);
    const double_v_t d2Psi_dIbar4Ibar4 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar4Ibar4Expression, parameterVector);
    const double_v_t d2Psi_dIbar1Ibar5 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar1Ibar5Expression, parameterVector);
    const double_v_t d2Psi_dIbar2Ibar5 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar2Ibar5Expression, parameterVector);
    const double_v_t d2Psi_dIbar5Ibar5 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar5Ibar5Expression, parameterVector);
    const double_v_t d2Psi_dIbar4Ibar5 = ExpressionHelper<double_v_t>::apply(d2Psi_dIbar4Ibar5Expression, parameterVector);

    // terms for 4th invariant
    decoupledFormFactor5 = 4*(d2Psi_dIbar1Ibar4 + Ibar1 * d2Psi_dIbar2Ibar4);
    decoupledFormFactor6 = -4*d2Psi_dIbar2Ibar4;
    decoupledFormFactor7 = 4*d2Psi_dIbar4Ibar4;

    // terms for 5th invariant
    decoupledFormFactor8 = 4*(d2Psi_dIbar1Ibar5 + Ibar1 * d2Psi_dIbar2Ibar5);
    decoupledFormFactor9 = -4*d2Psi_dIbar2Ibar5;
    decoupledFormFactor10 = 4*d2Psi_dIbar5Ibar5;
    decoupledFormFactor11 = 4*d2Psi_dIbar4Ibar5;
    decoupledFormFactor12 = 4*dPsi_dIbar5;
  }

  // prepare factors for coupled form of strain energy function, dependent on I1,I2,I3
  const double_v_t I1 = invariants[0];
  const double_v_t I3 = invariants[2];

  auto dPsi_dI1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupled, Term::I1);
  auto dPsi_dI2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupled, Term::I2);
  auto dPsi_dI3Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupled, Term::I3);

  auto d2Psi_dI1I1Expression = SEMT::deriv_t(dPsi_dI1Expression, Term::I1);
  auto d2Psi_dI1I2Expression = SEMT::deriv_t(dPsi_dI1Expression, Term::I2);
  auto d2Psi_dI1I3Expression = SEMT::deriv_t(dPsi_dI1Expression, Term::I3);
  auto d2Psi_dI2I2Expression = SEMT::deriv_t(dPsi_dI2Expression, Term::I2);
  auto d2Psi_dI2I3Expression = SEMT::deriv_t(dPsi_dI2Expression, Term::I3);
  auto d2Psi_dI3I3Expression = SEMT::deriv_t(dPsi_dI3Expression, Term::I3);

  const double_v_t dPsi_dI2 = ExpressionHelper<double_v_t>::apply(dPsi_dI2Expression, parameterVector);
  const double_v_t dPsi_dI3 = ExpressionHelper<double_v_t>::apply(dPsi_dI3Expression, parameterVector);
  const double_v_t d2Psi_dI1I1 = ExpressionHelper<double_v_t>::apply(d2Psi_dI1I1Expression, parameterVector);
  const double_v_t d2Psi_dI1I2 = ExpressionHelper<double_v_t>::apply(d2Psi_dI1I2Expression, parameterVector);
  const double_v_t d2Psi_dI1I3 = ExpressionHelper<double_v_t>::apply(d2Psi_dI1I3Expression, parameterVector);
  const double_v_t d2Psi_dI2I2 = ExpressionHelper<double_v_t>::apply(d2Psi_dI2I2Expression, parameterVector);
  const double_v_t d2Psi_dI2I3 = ExpressionHelper<double_v_t>::apply(d2Psi_dI2I3Expression, parameterVector);
  const double_v_t d2Psi_dI3I3 = ExpressionHelper<double_v_t>::apply(d2Psi_dI3I3Expression, parameterVector);

  const double_v_t coupledFormFactor1 = 4*(d2Psi_dI1I1 + 2*I1*d2Psi_dI1I2 + dPsi_dI2 + MathUtility::sqr(I1)*d2Psi_dI2I2);
  const double_v_t coupledFormFactor2 = -4*(d2Psi_dI1I2 + I1*d2Psi_dI2I2);
  const double_v_t coupledFormFactor3 = 4*(I3*d2Psi_dI1I3 + I1*I3*d2Psi_dI2I3);
  const double_v_t coupledFormFactor4 = 4*d2Psi_dI2I2;
  const double_v_t coupledFormFactor5 = -4*I3*d2Psi_dI2I3;
  const double_v_t coupledFormFactor6 = 4*(I3*dPsi_dI3 + MathUtility::sqr(I3)*d2Psi_dI3I3);
  const double_v_t coupledFormFactor7 = -4*I3*dPsi_dI3;
  const double_v_t coupledFormFactor8 = -4*dPsi_dI2;


  double_v_t pTilde = pressure;    // for incompressible material pTilde = pressure
  // for compressible material, pTilde = p + J*dp/dJ

  // for compressible material, the stress, p, is given by a constitutive equation
  // for incompressible material, the stress, p, is an unknown that will be solved for
  if (!Term::isIncompressible)    // if compressible material
  {
    auto dPsi_dJExpression = SEMT::deriv_t(Term::strainEnergyDensityFunctionVolumetric, Term::J);
    auto dp_dJExpression   = SEMT::deriv_t(dPsi_dJExpression, Term::J);
    const double_v_t dPsi_dJ = ExpressionHelper<double_v_t>::apply(dPsi_dJExpression, parameterVector);
    const double_v_t dp_dJ = ExpressionHelper<double_v_t>::apply(dp_dJExpression, parameterVector);

    pressure = dPsi_dJ;
    pTilde = pressure + J * dp_dJ;
  }

  // compute dIbar5_dC
  Tensor2<3,double_v_t> dIbar5_dC;

  // row index
  for (int i=0; i<3; i++)
  {
    // compute a0 Cbar_i a0
    double_v_t a0Cbar_i = 0;
    for (int k=0; k<3; k++)
    {
      //a0Cbar_i += a0_k * Cbar_ki
      a0Cbar_i += fiberDirection[k] * factorJ23*rightCauchyGreen[i][k];
    }

    // column index
    for (int j=0; j<3; j++)
    {
      // compute Cbar a0
      double_v_t cBarA0_j = 0;
      for (int k=0; k<3; k++)
      {
        //cBarA0_j += Cbar_jk * a0_k
        cBarA0_j += factorJ23*rightCauchyGreen[k][j] * fiberDirection[k];
      }
      dIbar5_dC[j][i] = fiberDirection[i] * cBarA0_j + a0Cbar_i * fiberDirection[j];
    }
  }

  // formula for C_iso: Holzapfel "Nonlinear Solid Mechanics" p.255
  // formula for C_bar: Holzapfel "Nonlinear Solid Mechanics" p.262

  std::array<double_v_t,21> entriesFromC;     // the 21 distinct values of ∂^2Ψ(C)/∂C∂C, only if the formulation includes the Ψ(C) term (the number of 21 is because of minor and major symmetrics in CC)

  // if the formulation includes a Ψ(C) term, we need to add CC = 4 * ∂^2Ψ(C)/∂C∂C
  if (usesFormulationWithC)
  {
    // C = F^T F is symmetric, we only use the entries C11, C12, C13, C22, C23, C33
    auto dPsi_dC11Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C11);
    auto dPsi_dC12Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C12);
    auto dPsi_dC13Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C13);
    auto dPsi_dC22Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C22);
    auto dPsi_dC23Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C23);
    auto dPsi_dC33Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionCoupledDependentOnC, Term::C33);

    // form the symbolic derivatives
    auto d2Psi_dC11dC11Expression = SEMT::deriv_t(dPsi_dC11Expression, Term::C11);
    auto d2Psi_dC12dC11Expression = SEMT::deriv_t(dPsi_dC12Expression, Term::C11);
    auto d2Psi_dC13dC11Expression = SEMT::deriv_t(dPsi_dC13Expression, Term::C11);
    auto d2Psi_dC22dC11Expression = SEMT::deriv_t(dPsi_dC22Expression, Term::C11);
    auto d2Psi_dC23dC11Expression = SEMT::deriv_t(dPsi_dC23Expression, Term::C11);
    auto d2Psi_dC33dC11Expression = SEMT::deriv_t(dPsi_dC33Expression, Term::C11);

    auto d2Psi_dC12dC12Expression = SEMT::deriv_t(dPsi_dC12Expression, Term::C12);
    auto d2Psi_dC13dC12Expression = SEMT::deriv_t(dPsi_dC13Expression, Term::C12);
    auto d2Psi_dC22dC12Expression = SEMT::deriv_t(dPsi_dC22Expression, Term::C12);
    auto d2Psi_dC23dC12Expression = SEMT::deriv_t(dPsi_dC23Expression, Term::C12);
    auto d2Psi_dC33dC12Expression = SEMT::deriv_t(dPsi_dC33Expression, Term::C12);

    auto d2Psi_dC13dC13Expression = SEMT::deriv_t(dPsi_dC13Expression, Term::C13);
    auto d2Psi_dC22dC13Expression = SEMT::deriv_t(dPsi_dC22Expression, Term::C13);
    auto d2Psi_dC23dC13Expression = SEMT::deriv_t(dPsi_dC23Expression, Term::C13);
    auto d2Psi_dC33dC13Expression = SEMT::deriv_t(dPsi_dC33Expression, Term::C13);

    auto d2Psi_dC22dC22Expression = SEMT::deriv_t(dPsi_dC22Expression, Term::C22);
    auto d2Psi_dC23dC22Expression = SEMT::deriv_t(dPsi_dC23Expression, Term::C22);
    auto d2Psi_dC33dC22Expression = SEMT::deriv_t(dPsi_dC33Expression, Term::C22);

    auto d2Psi_dC23dC23Expression = SEMT::deriv_t(dPsi_dC23Expression, Term::C23);
    auto d2Psi_dC33dC23Expression = SEMT::deriv_t(dPsi_dC33Expression, Term::C23);

    auto d2Psi_dC33dC33Expression = SEMT::deriv_t(dPsi_dC33Expression, Term::C33);

    // evaluate the 2nd derivatives of Ψ(C) at the given values of C
    const double_v_t d2Psi_dC11dC11 = ExpressionHelper<double_v_t>::apply(d2Psi_dC11dC11Expression, parameterVector);
    const double_v_t d2Psi_dC12dC11 = ExpressionHelper<double_v_t>::apply(d2Psi_dC12dC11Expression, parameterVector);
    const double_v_t d2Psi_dC13dC11 = ExpressionHelper<double_v_t>::apply(d2Psi_dC13dC11Expression, parameterVector);
    const double_v_t d2Psi_dC22dC11 = ExpressionHelper<double_v_t>::apply(d2Psi_dC22dC11Expression, parameterVector);
    const double_v_t d2Psi_dC23dC11 = ExpressionHelper<double_v_t>::apply(d2Psi_dC23dC11Expression, parameterVector);
    const double_v_t d2Psi_dC33dC11 = ExpressionHelper<double_v_t>::apply(d2Psi_dC33dC11Expression, parameterVector);

    const double_v_t d2Psi_dC12dC12 = ExpressionHelper<double_v_t>::apply(d2Psi_dC12dC12Expression, parameterVector);
    const double_v_t d2Psi_dC13dC12 = ExpressionHelper<double_v_t>::apply(d2Psi_dC13dC12Expression, parameterVector);
    const double_v_t d2Psi_dC22dC12 = ExpressionHelper<double_v_t>::apply(d2Psi_dC22dC12Expression, parameterVector);
    const double_v_t d2Psi_dC23dC12 = ExpressionHelper<double_v_t>::apply(d2Psi_dC23dC12Expression, parameterVector);
    const double_v_t d2Psi_dC33dC12 = ExpressionHelper<double_v_t>::apply(d2Psi_dC33dC12Expression, parameterVector);

    const double_v_t d2Psi_dC13dC13 = ExpressionHelper<double_v_t>::apply(d2Psi_dC13dC13Expression, parameterVector);
    const double_v_t d2Psi_dC22dC13 = ExpressionHelper<double_v_t>::apply(d2Psi_dC22dC13Expression, parameterVector);
    const double_v_t d2Psi_dC23dC13 = ExpressionHelper<double_v_t>::apply(d2Psi_dC23dC13Expression, parameterVector);
    const double_v_t d2Psi_dC33dC13 = ExpressionHelper<double_v_t>::apply(d2Psi_dC33dC13Expression, parameterVector);

    const double_v_t d2Psi_dC22dC22 = ExpressionHelper<double_v_t>::apply(d2Psi_dC22dC22Expression, parameterVector);
    const double_v_t d2Psi_dC23dC22 = ExpressionHelper<double_v_t>::apply(d2Psi_dC23dC22Expression, parameterVector);
    const double_v_t d2Psi_dC33dC22 = ExpressionHelper<double_v_t>::apply(d2Psi_dC33dC22Expression, parameterVector);

    const double_v_t d2Psi_dC23dC23 = ExpressionHelper<double_v_t>::apply(d2Psi_dC23dC23Expression, parameterVector);
    const double_v_t d2Psi_dC33dC23 = ExpressionHelper<double_v_t>::apply(d2Psi_dC33dC23Expression, parameterVector);

    const double_v_t d2Psi_dC33dC33 = ExpressionHelper<double_v_t>::apply(d2Psi_dC33dC33Expression, parameterVector);

    // distinct entries of C: {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
    // {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
    // {2,2,2,2}

    entriesFromC =
    {
      d2Psi_dC11dC11, d2Psi_dC12dC11, d2Psi_dC13dC11, d2Psi_dC22dC11, d2Psi_dC23dC11, d2Psi_dC33dC11,
      d2Psi_dC12dC12, d2Psi_dC13dC12, d2Psi_dC22dC12, d2Psi_dC23dC12, d2Psi_dC33dC12,
      d2Psi_dC13dC13, d2Psi_dC22dC13, d2Psi_dC23dC13, d2Psi_dC33dC13,
      d2Psi_dC22dC22, d2Psi_dC23dC22, d2Psi_dC33dC22,
      d2Psi_dC23dC23, d2Psi_dC33dC23,
      d2Psi_dC33dC33,
    };
  }

  // output for debugging
  if (false)
  {
    LOG(DEBUG) << "elasticity tensor, Ψ: " << Term::strainEnergyDensityFunctionIsochoric;
    LOG(DEBUG) << "∂Ψ/∂Ibar1: " << dPsi_dIbar1Expression;
    LOG(DEBUG) << "∂Ψ/∂Ibar2: " << dPsi_dIbar2Expression << " = " << dPsi_dIbar2;
    LOG(DEBUG) << "∂2Ψ/(∂Ibar1 ∂Ibar1): " << d2Psi_dIbar1Ibar1Expression << " = " << d2Psi_dIbar1Ibar1;
    LOG(DEBUG) << "∂2Ψ/(∂Ibar1 ∂Ibar1): " << d2Psi_dIbar1Ibar2Expression << " = " << d2Psi_dIbar1Ibar2;
    LOG(DEBUG) << "∂2Ψ/(∂Ibar1 ∂Ibar1): " << d2Psi_dIbar2Ibar2Expression << " = " << d2Psi_dIbar2Ibar2;
    LOG(DEBUG) << "decoupledFormFactor1: " << decoupledFormFactor1;
    LOG(DEBUG) << "decoupledFormFactor2: " << decoupledFormFactor2;
    LOG(DEBUG) << "decoupledFormFactor3: " << decoupledFormFactor3;
    LOG(DEBUG) << "decoupledFormFactor4: " << decoupledFormFactor4;
  }

  // distinct entries, only those have to be computed as the rest is symmetric
  const int indices[21][4] = {
    {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
    {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
    {2,2,2,2}
  };

  // to compute all entries and verify the symmetry, use the following code (needs further adjustments at  the loop boundaries: "entryNo<21" -> "entryNo<81")
#if 0
  int indices[81][4];

  int entryNo = 0;
  for (int a = 0; a < D; a++)
  {
    for (int b = 0; b < D; b++)
    {
      for (int c = 0; c < D; c++)
      {
        for (int d = 0; d < D; d++, entryNo++)
        {
          indices[entryNo][0] = a;
          indices[entryNo][1] = b;
          indices[entryNo][2] = c;
          indices[entryNo][3] = d;
        }
      }
    }
  }
  assert(entryNo == 81);
#endif

  // loop over distinct entries in elasticity tensor
  for (int entryNo = 0; entryNo < 21; entryNo++)
  {
    // get indices of current entry
    const int a = indices[entryNo][0];
    const int b = indices[entryNo][1];
    const int c = indices[entryNo][2];
    const int d = indices[entryNo][3];

    const double_v_t cInv_ab = inverseRightCauchyGreen[b][a];
    const double_v_t cInv_cd = inverseRightCauchyGreen[d][c];

    // compute C_vol
    const double_v_t cInvDotCInv = 0.5 * (inverseRightCauchyGreen[c][a] * inverseRightCauchyGreen[d][b] + inverseRightCauchyGreen[d][a] * inverseRightCauchyGreen[c][b]);

    const double_v_t cInvDyadCInv = inverseRightCauchyGreen[b][a] * inverseRightCauchyGreen[d][c];
    const double_v_t jp = deformationGradientDeterminant * pressure;
    const double_v_t jpTilde = deformationGradientDeterminant * pTilde;

    const double_v_t cVol = jpTilde * cInvDyadCInv - 2 * jp * cInvDotCInv;

    // compute C_iso
    //                          ab  ef     gh   cd
    // compute contribution from  P : Cbar : P^T
    double_v_t pCbarPT_abcd = 0.0;

    for (int g = 0; g < D; g++)
    {
      for (int h = 0; h < D; h++)
      {
        int delta_gh = (g == h? 1 : 0);

        double_v_t pCbar_abgh = 0;
        const double_v_t c_gh = rightCauchyGreen[h][g];

        for (int e = 0; e < D; e++)
        {
          int delta_eh = (e == h? 1 : 0);
          int delta_eg = (e == g? 1 : 0);

          for (int f = 0; f < D; f++)
          {
            int delta_fg = (f == g? 1 : 0);
            int delta_fh = (f == h? 1 : 0);
            int delta_ef = (e == f? 1 : 0);
            const double_v_t c_ef = rightCauchyGreen[f][e];

            int iI_efgh = delta_eg * delta_fh;
            int iIbar_efgh = delta_eh * delta_fg;

            // symmetric version
            double_v_t sS_efgh = 0.5 * (iI_efgh + iIbar_efgh);

            const double_v_t summand1 = decoupledFormFactor1 * delta_ef * delta_gh;
            const double_v_t summand2 = decoupledFormFactor2 * (delta_ef * factorJ23*c_gh + factorJ23*c_ef * delta_gh);
            const double_v_t summand3 = decoupledFormFactor3 * factorJ23*c_ef * factorJ23*c_gh;
            const double_v_t summand4 = decoupledFormFactor4 * sS_efgh;

            double_v_t sum = summand1 + summand2 + summand3 + summand4;

            // addd terms for the 4th and 5th invariants
            if (Term::usesFiberDirection)
            {
              // terms for 4th and 5th invariant
              const double_v_t summand5 = decoupledFormFactor5 * (delta_ef       * fiberDirection[g] * fiberDirection[h] + fiberDirection[e] * fiberDirection[f] * delta_gh);
              const double_v_t summand6 = decoupledFormFactor6 * (factorJ23*c_ef * fiberDirection[g] * fiberDirection[h] + fiberDirection[e] * fiberDirection[f] * factorJ23*c_gh);
              const double_v_t summand7 = decoupledFormFactor7 * fiberDirection[e] * fiberDirection[f] * fiberDirection[g] * fiberDirection[h];

              // terms for 5th invariant
              // indices efgh
              const double_v_t dIbar5_dC_ef = dIbar5_dC[f][e];
              const double_v_t dIbar5_dC_gh = dIbar5_dC[h][g];

              const double_v_t summand8  = decoupledFormFactor8 * (delta_ef * dIbar5_dC_gh + dIbar5_dC_ef * delta_gh);
              const double_v_t summand9  = decoupledFormFactor9 * (factorJ23*c_ef * dIbar5_dC_gh + dIbar5_dC_ef * factorJ23*c_gh);
              const double_v_t summand10 = decoupledFormFactor10 * (dIbar5_dC_ef * dIbar5_dC_gh);
              const double_v_t summand11 = decoupledFormFactor11 * (fiberDirection[e] * fiberDirection[f] * dIbar5_dC_gh + dIbar5_dC_ef * fiberDirection[g] * fiberDirection[h]);


              // ∂^2Ibar_5 / ∂Cbar_ef ∂Cbar_gh = 1/2 (δ_fh*a_e*a_g + δ_eh*a_f*a_g + δ_fg*a_e*a_h + δ_eg*a_f*a_h)
              double_v_t d2Ibar5_dCdC_efgh =
                delta_fh * fiberDirection[e] * fiberDirection[g]
                  + delta_eh * fiberDirection[f] * fiberDirection[g]
                  + delta_fg * fiberDirection[e] * fiberDirection[h]
                  + delta_eg * fiberDirection[f] * fiberDirection[h];

              const double_v_t summand12 = decoupledFormFactor12 * d2Ibar5_dCdC_efgh;

              sum += summand5 + summand6 + summand7;
              sum += summand8 + summand9 + summand10 + summand11 + summand12;
            }

            const double_v_t ccBar_efgh = factorJ43 * sum;

            //LOG(DEBUG) << "     CCBar_" << e << f << g << h << ": " << factorJ43 << "*(" << summand1 << "+" << summand2 << "+" << summand3 << "+" << summand4 << "), summand4 = " << factor4 << "*" << sS_efgh;

            fictitiousElasticityTensor[h][g][f][e] = ccBar_efgh;

            int delta_ae = (a == e? 1 : 0);
            int delta_bf = (b == f? 1 : 0);

            const double_v_t iI_abef = delta_ae * delta_bf;
            const double_v_t p_abef = iI_abef - 1./D*cInv_ab*c_ef;

            pCbar_abgh += p_abef * ccBar_efgh;

          }  // f
        }  // e

        int delta_cg = (c == g? 1 : 0);
        int delta_dh = (d == h? 1 : 0);
        const double_v_t iI_cdgh = delta_cg * delta_dh;

        const double_v_t p_cdgh = iI_cdgh - 1./D*cInv_cd*c_gh;     // (P^T)_ghcd = P_cdgh

        const double_v_t pT_ghcd = p_cdgh;

        pCbarPT_abcd += pCbar_abgh * pT_ghcd;

      }  // h
    }  // g

    // compute contribution from  2/3*J^{-2/3}*Sbar : C P_tilde
    double_v_t sBarC_abcd = 0.;
    for (int g = 0; g < D; g++)
    {
      for (int h = 0; h < D; h++)
      {
        sBarC_abcd += fictitiousPK2Stress[h][g] * rightCauchyGreen[h][g];
      }  // h
    }  // g

    const double_v_t pTilde_abcd = cInvDotCInv - 1./3 * cInvDyadCInv;

    const double_v_t sBarCP_abcd = 2./3. * factorJ23 * sBarC_abcd * pTilde_abcd;

    // compute contribution from -2./3*(CInv dyad Siso + Siso dyad CInv)
    const double_v_t cInvSiso = -2./3 * (cInv_ab * pk2StressIsochoric[d][c] + pk2StressIsochoric[b][a] * cInv_cd);

    // compute C_iso
    const double_v_t cIso = pCbarPT_abcd + sBarCP_abcd + cInvSiso;

    //LOG(DEBUG) << a<<b<<c<<d<<": compute CCIso: " << cIso << ": " << pCbarPT_abcd << "," << sBarCP_abcd << "," << cInvSiso;

#if 0
    // store entry in CCIso, this is only returned from this method for debugging purposes, only for use in materialTesting
    // therefore it can be omitted here
    elasticityTensorIso[d][c][b][a] = cIso;
    elasticityTensorIso[c][d][b][a] = cIso;
    elasticityTensorIso[d][c][a][b] = cIso;
    elasticityTensorIso[c][d][a][b] = cIso;
    elasticityTensorIso[b][a][d][c] = cIso;
    elasticityTensorIso[a][b][d][c] = cIso;
    elasticityTensorIso[b][a][c][d] = cIso;
    elasticityTensorIso[a][b][c][d] = cIso;
#endif
    // compute C
    double_v_t c_abcd = cVol + cIso;

    // debugging
    //const double_v_t c_abcd = cIso;

    // if the formulation includes a Ψ(C) term, we need to add CC = 4 * ∂^2Ψ(C)/∂C∂C
    if (usesFormulationWithC)
    {
      // add contribution to CC = 4 * ∂^2 Ψ(C)/∂C/∂C
      c_abcd += 4 * entriesFromC[entryNo];
    }

    // coupled form of strain energy
    // compute CC += δ1*I dyad I + ... + δ8*SS   (Holzapfel, p.272)

    int delta_ab = (a == b? 1 : 0);
    int delta_ac = (a == c? 1 : 0);
    int delta_ad = (a == d? 1 : 0);
    int delta_bc = (b == c? 1 : 0);
    int delta_bd = (b == d? 1 : 0);
    int delta_cd = (c == d? 1 : 0);
    const double_v_t c_ab = rightCauchyGreen[b][a];
    const double_v_t c_cd = rightCauchyGreen[d][c];

    int iI_abcd = delta_ac * delta_bd;
    int iIbar_abcd = delta_ad * delta_bc;

    // symmetric version
    const double_v_t sS_abcd = 0.5 * (iI_abcd + iIbar_abcd);

    const double_v_t summand1 = coupledFormFactor1 * delta_ab * delta_cd;
    const double_v_t summand2 = coupledFormFactor2 * (delta_ab * c_cd + c_ab * delta_cd);
    const double_v_t summand3 = coupledFormFactor3 * (delta_ab * cInv_cd + cInv_ab * delta_cd);
    const double_v_t summand4 = coupledFormFactor4 * c_ab * c_cd;
    const double_v_t summand5 = coupledFormFactor5 * (c_ab * cInv_cd + cInv_ab * c_cd);
    const double_v_t summand6 = coupledFormFactor6 * cInv_ab * cInv_cd;
    const double_v_t summand7 = coupledFormFactor7 * cInvDotCInv;
    const double_v_t summand8 = coupledFormFactor8 * sS_abcd;

    double_v_t sum = summand1 + summand2 + summand3 + summand4 + summand5 + summand6 + summand7 + summand8;

    // add contribution from coupled form of (compressible) strain energy function
    c_abcd += sum;

    // store entry, store all symmetric entries at once, because tensor has major and minor symmetries (CC_abcd = CC_cdab and CC_abcd = CC_bacd)
    elasticityTensor[d][c][b][a] = c_abcd;
    elasticityTensor[c][d][b][a] = c_abcd;
    elasticityTensor[d][c][a][b] = c_abcd;
    elasticityTensor[c][d][a][b] = c_abcd;
    elasticityTensor[b][a][d][c] = c_abcd;
    elasticityTensor[a][b][d][c] = c_abcd;
    elasticityTensor[b][a][c][d] = c_abcd;
    elasticityTensor[a][b][c][d] = c_abcd;
  }

  // verify major and minor symmetries of elasticity tensor
  if (VLOG_IS_ON(1))
  {
  //LOG_N_TIMES(2,DEBUG) << "elasticityTensor: " << elasticityTensor;

    const double_v_t errorTolerance = 1e-12;
    std::string name("C");
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            if (Vc::any_of(MathUtility::abs(elasticityTensor[a][b][c][d] - elasticityTensor[b][a][c][d]) > errorTolerance))
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" <<b<< "][" <<a<< "][" << c<< "][" << d<< "] (" <<elasticityTensor[b][a][c][d]<< " != " <<elasticityTensor[a][b][c][d]<< ") - minor symmetry violated" << std::endl;
            }
            if (Vc::any_of(MathUtility::abs(elasticityTensor[a][b][c][d] - elasticityTensor[a][b][d][c]) > errorTolerance))
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" <<a<< "][" <<b<< "][" << d<< "][" << c<< "] (" <<elasticityTensor[a][b][d][c]<< " != " <<elasticityTensor[a][b][c][d]<< ") - minor symmetry violated" << std::endl;
            }
            if (Vc::any_of(MathUtility::abs(elasticityTensor[a][b][c][d] - elasticityTensor[c][d][a][b]) > errorTolerance))
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" << c<< "][" << d<< "][" <<a<< "][" <<b<< "] (" <<elasticityTensor[c][d][a][b]<< " != " <<elasticityTensor[a][b][c][d]<< ") - major symmetry violated" << std::endl;
            }
          }
        }
      }
    }
    LOG_N_TIMES(2,DEBUG) << "elasticity tensor checked for minor symmetries";
    //LOG(DEBUG) << "elasticity tensor checked for minor symmetries";
  }
}

template<typename Term,typename MeshType, int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computePSbar(const Tensor2<3,double_v_t> &fictitiousPK2Stress, const Tensor2<3,double_v_t> &rightCauchyGreen)
{
  // only needed for debugging in materialTesting
  double_v_t determinant;
  Tensor2<3,double_v_t> inverseRightCauchyGreen = MathUtility::computeInverse<3>(rightCauchyGreen, determinant);

  //ab cd
  //  P : Sbar
  Tensor2<3,double_v_t> PSbar;
  for (int a = 0; a < 3; a++)
  {
    for (int b = 0; b < 3; b++)
    {
      double_v_t cInv_ab = inverseRightCauchyGreen[b][a];

      double_v_t psbar_ab = 0;
      for (int c = 0; c < 3; c++)
      {
        int delta_ac = (a == c? 1 : 0);
        for (int d = 0; d < 3; d++)
        {
          int delta_bd = (b == d? 1 : 0);
          const double_v_t ii_abcd = delta_ac * delta_bd;

          double_v_t c_cd = rightCauchyGreen[d][c];

          double_v_t p_abcd = ii_abcd - 1./3*cInv_ab * c_cd;

          //LOG(DEBUG) << "  .PP_" << a << b << c << d << ": " << p_abcd << "=" << ii_abcd << "-1/3*" << cInv_ab * c_cd << " (" << cInv_ab << "," << c_cd << "), Ii: " << ii_abcd;

          const double_v_t sBar_cd = fictitiousPK2Stress[d][c];
          psbar_ab += p_abcd * sBar_cd;
        }
      }


      //LOG(DEBUG) << " PSBar_" << a << b << ": " << psbar_ab;
      PSbar[b][a] = psbar_ab;
    }
  }
  return PSbar;
}

} // namespace

