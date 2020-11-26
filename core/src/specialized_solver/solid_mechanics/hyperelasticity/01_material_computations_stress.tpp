#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include <Vc/Vc>

#include "equation/mooney_rivlin_incompressible.h"
#include "utility/math_utility.h"
#include "utility/vector_operators.h"

namespace SpatialDiscretization
{

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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

  // Ibar1, Ibar2, Ibar4, Ibar5, J, I1, I2, I3, C11, C12, C13, C22, C23, C33, a1, a2, a3
  std::vector<double_v_t> parameterVector = {
    reducedInvariants[0], reducedInvariants[1], reducedInvariants[3], reducedInvariants[4],  // Ibar1, Ibar2, Ibar4, Ibar5
    deformationGradientDeterminant,                                                          // J
    invariants[0], invariants[1], invariants[2],                                             // I1, I2, I3
    rightCauchyGreen[0][0], rightCauchyGreen[1][0], rightCauchyGreen[2][0],                  // C11, C12, C13
    rightCauchyGreen[1][1], rightCauchyGreen[2][1], rightCauchyGreen[2][2],                  // C22, C23, C33
    fiberDirection[0], fiberDirection[1], fiberDirection[2]                                  // a1, a2, a3
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

  // the following line gives linker errors in debug target
  //VLOG(2) << "coupled term: " << Term::strainEnergyDensityFunctionCoupled;
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
    //VLOG(1) << "J = " << Term::J << " = " << J << ", Psi_vol: " << Term::strainEnergyDensityFunctionVolumetric << ", dΨ/dJ: " << dPsi_dJExpression << " = " << dPsi_dJ;
    pressure = dPsi_dJ;
  }

  double factor23 = -2./3;
  double_v_t factorJ23 = MathUtility::pow(J, factor23);

  // If the jacobian is negative, this is a non-physical state. A warning has already been printed. Now to continue computation, set J^-2/3 to 1 (would be nan otherwise).
  if (Vc::any_of(J < 0))
    factorJ23 = 1;

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

#ifndef NDEBUG
  VLOG(1) << "fictitiousPK2Stress: " << fictitiousPK2Stress << ", inverseRightCauchyGreen: " << inverseRightCauchyGreen << ", J: " << J << ", p: " << pressure;
  VLOG(1) << "decoupledFormFactors: " << decoupledFormFactor1 << ", " << decoupledFormFactor2 << ", " << decoupledFormFactor4 << ", " << decoupledFormFactor5;
  VLOG(1) << "fiberDirection: " << fiberDirection << ", C: " << rightCauchyGreen;
#endif

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
      //VLOG(2) << "   sVol: " << sVol << ", sIso: " << sIso;

      // total stress is sum of volumetric and isochoric part, sVol = J*p*C^{-1}_AB, sIso = j^{-2/3}*(Ii-1/3*Cc)*Sbar
      pK2Stress[j][i] = sVol + sIso;

      // debugging
      //pK2Stress[j][i] = sIso;

      VLOG(2) << "set pk2Stress_" << i << j << " = " << pK2Stress[j][i];
      //LOG(DEBUG) << "set pk2Stress_" << i << j << " = " << pK2Stress[j][i] << " = " << sVol << " + " << sIso;

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
  //LOG(DEBUG) << "   δ1: " << coupledFormFactor1 << ", δ2: " << coupledFormFactor2 << ", δ3: " << coupledFormFactor3 << ", after adding coupled form terms: " << pK2Stress;

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
    //LOG(DEBUG) << "formulation with C, add S=(" << dPsi_dC11 << "," << 2*dPsi_dC12 << "," << 2*dPsi_dC13 << "," << 2*dPsi_dC22 << "," << 2*dPsi_dC23 << "," << 2*dPsi_dC33 << "), parameterVector: " << parameterVector;
    //LOG(DEBUG) << "expressions: (" << dPsi_dC11Expression << "\n2: " << dPsi_dC12Expression << "\n3: " << dPsi_dC13Expression << "\n4: " << dPsi_dC22Expression << "\n5: " << dPsi_dC23Expression << "\n6: " << dPsi_dC33Expression << ")";

    // add contribution to S = 2 * ∂Ψ(C)/∂C
    pK2Stress[0][0] += 2*dPsi_dC11;   // S11          // S11
    pK2Stress[1][0] += 2*dPsi_dC12;   // S12          // S12
    pK2Stress[2][0] += 2*dPsi_dC13;   // S13          // S13
    pK2Stress[0][1] += 2*dPsi_dC12;   // S21 = S12    // S21
    pK2Stress[1][1] += 2*dPsi_dC22;   // S22          // S22
    pK2Stress[2][1] += 2*dPsi_dC23;   // S23          // S23
    pK2Stress[0][2] += 2*dPsi_dC13;   // S31 = S13    // S31
    pK2Stress[1][2] += 2*dPsi_dC23;   // S32 = S23    // S32
    pK2Stress[2][2] += 2*dPsi_dC33;   // S33          // S33
  }
  /*static int cntr = -1;
  cntr++;
  LOG(DEBUG) << "pK2Stress: " << pK2Stress << " counter: " << cntr;
  */

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
    if (MathUtility::containsNanOrInf(fictitiousPK2Stress))
    {
      LOG(ERROR) << "fictitiousPK2Stress contains nan: " << fictitiousPK2Stress << ", J=" << J;
    }
    LOG(FATAL) << "PK2stress contains nan: " << pK2Stress;
  }
#endif

  return pK2Stress;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
computePK2StressField()
{
  //LOG(TRACE) << "computePK2StressField";

  //this->data_.pK2Stress()->startGhostManipulation();
  this->data_.pK2Stress()->zeroGhostBuffer();
  this->data_.materialTraction()->zeroGhostBuffer();
  this->data_.materialTraction()->zeroEntries();

  this->data_.traction()->zeroGhostBuffer();
  this->data_.traction()->zeroEntries();

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
    double_v_t approximateMeshWidth = MathUtility::computeApproximateMeshWidth<double_v_t,nDisplacementsDofsPerElement>(geometryReferenceValues);

    // get displacements field values for current element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> displacementsValues;
    this->data_.displacements()->getElementValues(elementNoLocalv, displacementsValues);

    // get velocity field values for current element
    std::array<Vec3_v_t,nDisplacementsDofsPerElement> velocitiesValues;
    this->data_.velocities()->getElementValues(elementNoLocalv, velocitiesValues);

    // get pressure field values for current element
    std::array<double_v_t,nPressureDofsPerElement> pressureValuesCurrentElement;
    this->data_.pressure()->getElementValues(elementNoLocalv, pressureValuesCurrentElement);

    // get direction values for current element
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
        inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, approximateMeshWidth, jacobianDeterminant);

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
          // row-major, F_ij
          deformationGradientValues[i*3+j] = deformationGradient[j][i];
          deformationGradientTimeDerivativeValues[i*3+j] = Fdot[j][i];
        }
      }
      this->data_.deformationGradient()->setValue(dofNoLocal, deformationGradientValues, INSERT_VALUES);
      this->data_.deformationGradientTimeDerivative()->setValue(dofNoLocal, deformationGradientTimeDerivativeValues, INSERT_VALUES);


      Tensor2_v_t<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F

      double_v_t rightCauchyGreenDeterminant;   // J^2
      Tensor2_v_t<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, approximateMeshWidth, rightCauchyGreenDeterminant);  // C^-1

      // fiber direction
      Vec3_v_t fiberDirection = displacementsFunctionSpace->template interpolateValueInElement<3>(elementalDirectionValues, xi);

      // fiberDirection is not automatically normalized because of the interpolation inside the element, normalize again
      if (Term::usesFiberDirection)
      {
        MathUtility::normalize<3>(fiberDirection);
      }

#ifndef NDEBUG
      if (Term::usesFiberDirection)
      {
        if (fabs(MathUtility::norm<3>(fiberDirection) - 1) > 1e-3)
          LOG(FATAL) << "fiberDirecton " << fiberDirection << " is not normalized (c)(norm: " << MathUtility::norm<3>(fiberDirection)
            << ", difference to 1: " << MathUtility::norm<3>(fiberDirection) - 1 << ") elementalDirectionValues:" << elementalDirectionValues;
      }
#endif

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
        LOG(FATAL) << "reducedInvariants contains nan: " << reducedInvariants
          << ", invariants: " << invariants << ", deformationGradient: " << deformationGradient
          << ", deformationGradientDeterminant: " << deformationGradientDeterminant
          << ", rightCauchyGreenDeterminant: " << rightCauchyGreenDeterminant << ", rightCauchyGreen: " << rightCauchyGreen;

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
      Mesh::face_t outwardsFace = Mesh::face_t::face2Minus;

      // get normal
      if (indexZ == 0 || indexZ == 1)
      {
        // bottom and center nodes of a quadratic element
        outwardsFace = Mesh::face_t::face2Minus;
      }
      else if (indexZ == 2)
      {
        // top node of a quadratic element
        outwardsFace = Mesh::face_t::face2Plus;
      }

      // bottom node, get normal in reference configuration
      Vec3_v_t materialNormal = displacementsFunctionSpace->getNormal(outwardsFace, geometryReferenceValues, xi);

      // get normal in current configuration
      //Vec3_v_t normal = displacementsFunctionSpace->getNormal(outwardsFace, elementNoLocalv, xi);

      // compute traction by Cauchy theorem T = S n0
      Vec3_v_t materialTraction = pK2Stress * materialNormal;

      //LOG(INFO) << elementNoLocal << "." << elementalNodeNo << ": traction: " << traction;

      // set value in material traction
      this->data_.materialTraction()->setValue(dofNoLocal, materialTraction, INSERT_VALUES);

      // S n0 = F^-1 t  with S=symmetric PK2 stress, n0=normal in reference configuration, t=traction in current configuration
      // => t = F S n0

      // compute traction in current configuration
      Vec3_v_t traction = deformationGradient * materialTraction;

      // set value in traction
      this->data_.traction()->setValue(dofNoLocal, traction, INSERT_VALUES);
    }
  }
  this->data_.pK2Stress()->zeroGhostBuffer();
  this->data_.pK2Stress()->finishGhostManipulation();
  this->data_.pK2Stress()->startGhostManipulation();

  this->data_.materialTraction()->zeroGhostBuffer();
  this->data_.materialTraction()->finishGhostManipulation();
  this->data_.materialTraction()->startGhostManipulation();

  this->data_.traction()->zeroGhostBuffer();
  this->data_.traction()->finishGhostManipulation();
  this->data_.traction()->startGhostManipulation();

  this->data_.deformationGradient()->zeroGhostBuffer();
  this->data_.deformationGradient()->finishGhostManipulation();
  this->data_.deformationGradient()->startGhostManipulation();

  this->data_.deformationGradientTimeDerivative()->zeroGhostBuffer();
  this->data_.deformationGradientTimeDerivative()->finishGhostManipulation();
  this->data_.deformationGradientTimeDerivative()->startGhostManipulation();
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
computeBearingForceAndMoment(const std::vector<std::tuple<element_no_t,bool>> &elements,
                             Vec3 &bearingForceBottom, Vec3 &bearingMomentBottom, Vec3 &bearingForceTop, Vec3 &bearingMomentTop)
{
  // set result variables to zero
  bearingForceBottom = Vec3{0};
  bearingMomentBottom = Vec3{0};
  bearingForceTop = Vec3{0};
  bearingMomentTop = Vec3{0};

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D-1,Quadrature::Gauss<2>> QuadratureSurface;
  typedef Vec3 EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureSurface::numberEvaluations()
          > EvaluationsArrayType;

  // setup arrays used for integration
  std::array<std::array<double,D-1>, QuadratureSurface::numberEvaluations()> samplingPoints = QuadratureSurface::samplingPoints();
  EvaluationsArrayType evaluationsArrayForce{};
  EvaluationsArrayType evaluationsArrayMoment{};

  for (int elementNoIndex = 0; elementNoIndex < elements.size(); elementNoIndex++)
  {
    element_no_t elementNoLocal = std::get<0>(elements[elementNoIndex]);
    bool isTop = std::get<1>(elements[elementNoIndex]);

    // get geometry field of reference configuration
    std::array<Vec3,nDisplacementsDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference()->getElementValues(elementNoLocal, geometryReferenceValues);

    // get traction t values in current element
    std::array<Vec3,nDisplacementsDofsPerElement> tractionValues;
    this->data_.traction()->getElementValues(elementNoLocal, tractionValues);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D-1> xiSurface = samplingPoints[samplingPointIndex];
      double xi2 = 0;
      if (isTop)
        xi2 = 1;
      std::array<double,D> xi{xiSurface[0], xiSurface[1], xi2};

      // compute the 3xD jacobian of the parameter space to world space mapping
      std::array<Vec3,D> jacobianMaterial = this->displacementsFunctionSpace_->computeJacobian(geometryReferenceValues, xi);

      // compute the traction value at the current sampling point xi
      Vec3 traction = this->displacementsFunctionSpace_->interpolateValueInElement(tractionValues, xi);

      // compute the position in reference configuration of the current sampling point
      Vec3 point = this->displacementsFunctionSpace_->interpolateValueInElement(geometryReferenceValues, xi);

      double integrationFactor = MathUtility::computeIntegrationFactor(jacobianMaterial);
      evaluationsArrayForce[samplingPointIndex] = traction * integrationFactor;
      evaluationsArrayMoment[samplingPointIndex] = traction * point * integrationFactor;
    }

    // integrate all values in the element at once
    EvaluationsType integratedValuesForce = QuadratureSurface::computeIntegral(evaluationsArrayForce);
    EvaluationsType integratedValuesMoment = QuadratureSurface::computeIntegral(evaluationsArrayMoment);

    if (isTop)
    {
      bearingForceTop += integratedValuesForce;
      bearingMomentTop += integratedValuesMoment;
    }
    else
    {
      bearingForceBottom += integratedValuesForce;
      bearingMomentBottom += integratedValuesMoment;
    }
  }

  Vec3 bearingForceBottomLocal = bearingForceBottom;
  Vec3 bearingMomentBottomLocal = bearingMomentBottom;
  Vec3 bearingForceTopLocal = bearingForceTop;
  Vec3 bearingMomentTopLocal = bearingMomentTop;

  // reduce values
  MPI_Comm mpiCommunicator = this->displacementsFunctionSpace_->meshPartition()->rankSubset()->mpiCommunicator();
  MPIUtility::handleReturnValue(MPI_Allreduce(bearingForceBottomLocal.data(), bearingForceBottom.data(),
                                              3, MPI_DOUBLE, MPI_SUM, mpiCommunicator), "MPI_Allreduce");

  MPIUtility::handleReturnValue(MPI_Allreduce(bearingMomentBottomLocal.data(), bearingMomentBottom.data(),
                                              3, MPI_DOUBLE, MPI_SUM, mpiCommunicator), "MPI_Allreduce");

  MPIUtility::handleReturnValue(MPI_Allreduce(bearingForceTopLocal.data(), bearingForceTop.data(),
                                              3, MPI_DOUBLE, MPI_SUM, mpiCommunicator), "MPI_Allreduce");

  MPIUtility::handleReturnValue(MPI_Allreduce(bearingMomentTopLocal.data(), bearingMomentTop.data(),
                                              3, MPI_DOUBLE, MPI_SUM, mpiCommunicator), "MPI_Allreduce");


}

} // namespace

