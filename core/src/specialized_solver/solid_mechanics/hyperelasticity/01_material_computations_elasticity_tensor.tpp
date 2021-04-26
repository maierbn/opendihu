#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available

#include "equation/mooney_rivlin_incompressible.h"
#include "utility/math_utility.h"

namespace SpatialDiscretization
{

//! compute the material elasticity tensor
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
void HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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

  // if the alternative coupled form of the strain energy function, ψ(C,a), strainEnergyDensityFunctionCoupledDependentOnC, is considered
  const bool usesFormulationWithC = typeid(decltype(Term::strainEnergyDensityFunctionCoupledDependentOnC)) != typeid(decltype(INT(0)));

  // Ibar1, Ibar2, Ibar4, Ibar5, J, I1, I2, I3, C11, C12, C13, C22, C23, C33, a1, a2, a3
  std::vector<double_v_t> parameterVector = {
    reducedInvariants[0], reducedInvariants[1], reducedInvariants[3], reducedInvariants[4],  // Ibar1, Ibar2, Ibar4, Ibar5
    deformationGradientDeterminant,                                                          // J
    invariants[0], invariants[1], invariants[2],                                             // I1, I2, I3
    rightCauchyGreen[0][0], rightCauchyGreen[1][0], rightCauchyGreen[2][0],                  // C11, C12, C13
    rightCauchyGreen[1][1], rightCauchyGreen[2][1], rightCauchyGreen[2][2],                  // C22, C23, C33
    fiberDirection[0], fiberDirection[1], fiberDirection[2]                                  // a1, a2, a3
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

    // distinct entries of C: {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},
    // {0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1}, {2,2,0,1},
    // {0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},
    // {1,1,1,1},{1,2,1,1},{2,2,1,1},
    // {1,2,1,2},{2,2,1,2},
    // {2,2,2,2}

    entriesFromC =
    {
      d2Psi_dC11dC11, d2Psi_dC12dC11, d2Psi_dC13dC11, d2Psi_dC22dC11, d2Psi_dC23dC11, d2Psi_dC33dC11,
      d2Psi_dC12dC12, d2Psi_dC13dC12, d2Psi_dC22dC12, d2Psi_dC23dC12, d2Psi_dC33dC12,
      d2Psi_dC13dC13, d2Psi_dC22dC13, d2Psi_dC23dC13, d2Psi_dC33dC13,
      d2Psi_dC22dC22, d2Psi_dC23dC22, d2Psi_dC33dC22,
      d2Psi_dC23dC23, d2Psi_dC33dC23,
      d2Psi_dC33dC33
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

            // add terms for the 4th and 5th invariants
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
#ifndef HAVE_STDSIMD
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
#endif
}

} // namespace

