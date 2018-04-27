#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

#include "easylogging++.h"
#include "utility/math_utility.h"
#include "control/types.h"

#include <cmath>
#include <array>

namespace SpatialDiscretization
{

template<typename BasisOnMeshType, typename Term>
Tensor2<BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeDeformationGradientParameterSpace(const std::array<Vec3,BasisOnMeshType::nDofsPerElement()> &displacement,
                                         const std::array<double, BasisOnMeshType::dim()> xi)
{
 
  // compute the deformation gradient w.r.t parameter space, dx/dxi = x_i,j = d_ij + u_i,j
  // where j is dimensionColumn and i is component of the used Vec3's
  // The real deformation gradient would be dx/dX, but computed here is dx/dxi 
  Tensor2<BasisOnMeshType::dim()> deformationGradientParameterSpace;
 
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < BasisOnMeshType::dim(); dimensionColumn++)
  {
    Vec3 du_dxi({0});   // handle full-dimension vector of displacement (i.e. (x,y,z))
    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      double dphi_dxi = BasisOnMeshType::dphi_dxi(dofIndex, dimensionColumn, xi);
      du_dxi += dphi_dxi * displacement[dofIndex];   // vector-valued addition
    }
    
    deformationGradientParameterSpace[dimensionColumn] = du_dxi;
    
    // add Kronecker delta to obtain x_i,j = delta_ij + u_i,j
    deformationGradientParameterSpace[dimensionColumn][dimensionColumn] += 1;
  }
  
  return deformationGradientParameterSpace;
}

template<typename BasisOnMeshType, typename Term>
Tensor2<BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeDeformationGradient(const std::array<Vec3,BasisOnMeshType::nDofsPerElement()> &displacement,
                           const Tensor2<BasisOnMeshType::dim()> &inverseJacobianMaterial,
                           const std::array<double, BasisOnMeshType::dim()> xi
                          )
{
  // compute the deformation gradient x_i,j = d_ij + u_i,j
  // where j is dimensionColumn and i is component of the used Vec3's
  
 
  VLOG(3) << "compute deformation gradient Fij, displacements: " << displacement;
 
  const int D = BasisOnMeshType::dim();
  Tensor2<D> deformationGradient;
 
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  // easier understandable version, less optimized
#if 0  
  VLOG(3) << "A=1, a=2";
  
  for (int A = 0; A < BasisOnMeshType::dim(); A++) // dimensionColumn
  {
    for (int a = 0; a < BasisOnMeshType::dim(); a++)
    {
      VLOG(3) << " entry " << a << A;
      double sum = 0;
      double sum_dphi_dxil2 = 0;
      for (int M = 0; M < nDofsPerElement; M++)
      {
        double dPhi_dXA = 0;
        for (int l = 0; l < BasisOnMeshType::dim(); l++)
        {
          double dphi_dxil = BasisOnMeshType::dphi_dxi(M, l, xi);
          double dxil_dXA = inverseJacobianMaterial[A][l];     // inverseJacobianMaterial[A][l] = J_lA = dxi_l/dX_A
        
          //VLOG(3) << "  L=" << M << ", l=" << l << ", dphiL/dxil=" << dphi_dxil << ", dxil/dXA=" << dxil_dXA << "u^L_a=" << displacement[M][a];
        
          dPhi_dXA += dphi_dxil * dxil_dXA;
        }
        VLOG(3) << "  L=" << M << ", dPhi_dXA: " << dPhi_dXA << ", u: " << displacement[M][a];
        
        sum += dPhi_dXA * displacement[M][a];
        sum_dphi_dxil2 += BasisOnMeshType::dphi_dxi(M, 2, xi);
      }
      VLOG(3) << "    sum_dphi_dxil2: " << sum_dphi_dxil2;
      deformationGradient[A][a] = (a==A? 1: 0) + sum;
      VLOG(3) << " entry " << a << A << " = " << deformationGradient[A][a];
    }
  }
#endif  

  // loop over dimension, i.e. columns of deformation gradient, j
  for (int dimensionColumn = 0; dimensionColumn < BasisOnMeshType::dim(); dimensionColumn++)
  {
    // compute du_i/dX_columnIndex
    std::array<double,D> du_dX({0});   // handle full-dimension vector of displacement (i.e. (x,y,z))
    
    VLOG(3) << " j = " << dimensionColumn;
    
    // loop over dof no., M
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      VLOG(3) << "  M = " << dofIndex;
      // compute dphi_dofIndex/dX_dimensionColumn
      double dphi_dX = 0;
      for (int l = 0; l < BasisOnMeshType::dim(); l++)
      {
        VLOG(3) << "   l = " << l;
        double dphi_dxil = BasisOnMeshType::dphi_dxi(dofIndex, l, xi);
        double dxil_dX = inverseJacobianMaterial[dimensionColumn][l];     // inverseJacobianMaterial[j][l] = J_lj = dxi_l/dX_j
        
        VLOG(3) << "     dphi_dxil = " << dphi_dxil << ", dxil_dX = " << dxil_dX;
        dphi_dX += dphi_dxil * dxil_dX;
      }
      VLOG(3) << "   dphi_dX = " << dphi_dX;
      
      // multiply dphi/dxi with dxi/dX to obtain dphi/dX
      VLOG(3) << "   displ_M: " << displacement[dofIndex];
      du_dX += dphi_dX * MathUtility::transformToD<D,3>(displacement[dofIndex]);   // vector-valued addition
    }
    VLOG(3) << " du_dXj: " << du_dX;
    
    deformationGradient[dimensionColumn] = du_dX;
    
    // add Kronecker delta to obtain x_i,j = delta_ij + u_i,j
    deformationGradient[dimensionColumn][dimensionColumn] += 1;
  }
  
  VLOG(3) << "deformationGradient: " << deformationGradient;
  return deformationGradient;
}

template<typename BasisOnMeshType, typename Term>
Tensor2<BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeRightCauchyGreenTensor(const Tensor2<BasisOnMeshType::dim()> &deformationGradient)
{
  // compute C = F^T*F where F is the deformationGradient and C is the right Cauchy-Green Tensor
  // the quantities are 3x3 tensors for the 3D case
  Tensor2<BasisOnMeshType::dim()> rightCauchyGreenTensor({std::array<double,BasisOnMeshType::dim()>({0})});
  
  // loop over dimension, i.e. columns of right cauchy green tensor, i
  for (int dimensionColumn = 0; dimensionColumn < BasisOnMeshType::dim(); dimensionColumn++)
  {
    
    // loop over row of tensor, j
    for (int dimensionRow = 0; dimensionRow < BasisOnMeshType::dim(); dimensionRow++)
    {
      for (int k = 0; k < BasisOnMeshType::dim(); k++)
      {
        rightCauchyGreenTensor[dimensionColumn][dimensionRow] += 
          deformationGradient[dimensionColumn][k] * deformationGradient[dimensionRow][k];
      }
    }
  }
  
  return rightCauchyGreenTensor;
}

template<typename BasisOnMeshType, typename Term>
std::array<double,3> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeInvariants(const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen, const double rightCauchyGreenDeterminant)
{
  std::array<double,3> invariants;
  
  // I1 = tr(C)
  invariants[0] = 0.0;
  for (int i = 0; i < BasisOnMeshType::dim(); i++)
  {
    invariants[0] += rightCauchyGreen[i][i];
  }
  
  // tr(C^2)
  double traceC2 = 0;
  for (int i = 0; i < BasisOnMeshType::dim(); i++)
  {
    for (int j = 0; j < BasisOnMeshType::dim(); j++)
    {
      traceC2 += rightCauchyGreen[i][j] * rightCauchyGreen[j][i];
    }
  }
  
  // I2 = 1/2 * (tr(C)^2 - tr(C^2))
  invariants[1] = 0.5 * (MathUtility::sqr(invariants[0]) - traceC2);
  
  // I3 = det(C)
  invariants[2] = rightCauchyGreenDeterminant;
  
  return invariants;
}


template<typename BasisOnMeshType, typename Term>
std::array<double,2> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeReducedInvariants(const std::array<double,3> invariants, const double deformationGradientDeterminant)
{
  std::array<double,2> reducedInvariants;
  
  // 3D: Fbar = J^{1/3}*F such that det(Fbar) = 1
  // 2D: Fbar = J^{1/2}*F such that det(Fbar) = 1
  
  // 3D: Cbar = J^{2/3}*C such that det(Cbar) = 1
  // 2D: Cbar = J^{2/2}*C = J*C = J*F^T*F = J*J^{-1/2}*J^{-1/2}*Fbar^T*Fbar such that det(Cbar) = 1
  
  
  double factor23 = -2./3;
  double factor43 = -4./3;
  
  if (BasisOnMeshType::dim() == 2)
  {
    factor23 = -1./2;
    factor43 = -2./2;
  }
  
  reducedInvariants[0] = pow(deformationGradientDeterminant, factor23) * invariants[0];
  reducedInvariants[1] = pow(deformationGradientDeterminant, factor43) * invariants[1];
  
  return reducedInvariants;
}

template<typename BasisOnMeshType, typename Term>
double SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeArtificialPressure(const double deformationGradientDeterminant, double &artificialPressureTilde)
{
  //Psi_vol = kappa*1/2*(J-1)^2, dPsi_vol/dJ = kappa*(J-1)
 
  // compute the artifical pressure for penalty formulation, p = dPsi_vol/dJ 
  auto dPsidJExpression = SEMT::deriv_t(Term::strainEnergyDensityFunctionVolumetric, Term::J);
  auto dpdJExpression = SEMT::deriv_t(dPsidJExpression, Term::J);
  
  VLOG(2) << "computeArtificialPressure, psi=" << Term::strainEnergyDensityFunctionVolumetric << ", dPsi/dJ=" << dPsidJExpression << ", dp/dJ=" << dpdJExpression;
  
  std::vector<double> deformationGradientDeterminantVector(3,deformationGradientDeterminant);
  
  const double artificialPressure = dPsidJExpression.apply(deformationGradientDeterminantVector);
  const double dpdJ = dpdJExpression.apply(deformationGradientDeterminantVector);
  
  VLOG(2) << "  J = 1+" << (deformationGradientDeterminant-1) << " kappa=" << SEMT::Parameter<2>::get_value() << " => artificialPressure p = " << artificialPressure << ", dpdJ: " << dpdJ;
  
  // pTilde = p + J*dp/dJ
  artificialPressureTilde = artificialPressure + deformationGradientDeterminant*dpdJ;
  
  return artificialPressure;
}

template<typename BasisOnMeshType, typename Term>
Tensor2<BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computePK2Stress(const double pressure,                             //< [in] pressure value p
                 const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen,        //< [in] C
                 const Tensor2<BasisOnMeshType::dim()> &inverseRightCauchyGreen, //< [in] C^{-1}
                 const std::array<double,2> reducedInvariants,     //< [in] the reduced invariants Ibar_1, Ibar_2 
                 const double deformationGradientDeterminant,      //< [in] J = det(F)
                 Tensor2<BasisOnMeshType::dim()> &fictitiousPK2Stress,          //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                 Tensor2<BasisOnMeshType::dim()> &pk2StressIsochoric            //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                )
{
  // compute the PK2 stress tensor as S=2*dPsi/dC
  // for explanation see pdf document
  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);
  
  std::vector<double> reducedInvariantsVector(reducedInvariants.begin(), reducedInvariants.end());
  
  const double dPsi_dIbar1 = dPsi_dIbar1Expression.apply(reducedInvariantsVector);
  const double dPsi_dIbar2 = dPsi_dIbar2Expression.apply(reducedInvariantsVector);
  
  const double Ibar1 = reducedInvariants[0];
  const double J = deformationGradientDeterminant;
  
  double factor23 = -2./3;
  if (BasisOnMeshType::dim() == 2)
  {
    factor23 = -1./2;
  }
  
  double factor1 = 2*(dPsi_dIbar1 + Ibar1*dPsi_dIbar2);
  double factor2 = -2*dPsi_dIbar2;
  double factorJ23 = pow(J, factor23);
  
  Tensor2<BasisOnMeshType::dim()> pK2Stress;
  
  // compute fictitiousPK2Stress:
  // Sbar = factor1*I + factor2*Cbar, Cbar = J^{-2/3}*C
  
  // row index
  for (int i=0; i<BasisOnMeshType::dim(); i++)
  {
    // column index
    for (int j=0; j<BasisOnMeshType::dim(); j++)
    {
      int delta_ij = (i == j? 1 : 0);
      const double cBar = factorJ23 * rightCauchyGreen[j][i];
      fictitiousPK2Stress[j][i] = factor1 * delta_ij + factor2 * cBar;
      
      //if (i == j)
      //  LOG(DEBUG) << "  C_" << i << j << " = " << rightCauchyGreen[j][i] << ", factorJ23 = " << factorJ23 << ", SBar_" << i << j << ": " << fictitiousPK2Stress[j][i] 
      //    << " = " << factor1 << " * " << delta_ij << " + " << factor2 << "*" << cBar ;
    }
  }
  
  // Holzapfel p.234
  // S = S_vol + S_iso
  // S_vol = J*p*C^-1
  // S_iso = 2*dPsi_iso/dC = J^(-2/3) P : Sbar   (P = II - 1/3(C^{-1} dyad C): Holzapfel p.229)
  // pSbar = P : Sbar, P: 4th order tensor, Sbar: 2nd order tensor, note: A:B = A_ijkl*B_kl*e_i dyad e_j
  
  // compute S = S_vol + S_iso
  // row index
  for (int i=0; i<BasisOnMeshType::dim(); i++)
  {
    // column index
    for (int j=0; j<BasisOnMeshType::dim(); j++)
    {
      // volumetric stress
      const double sVol = J * pressure * inverseRightCauchyGreen[j][i];
      
      // compute P : Sbar
      double pSbar = 0;
      //double ccs = 0;  // (C^-1 dyad C)*Sbar
      //double cSbar = 0;  // C:Sbar
      
      // row index
      for (int k=0; k<BasisOnMeshType::dim(); k++)
      {
        const int delta_ik = (i == k? 1 : 0);
        
        // column index
        for (int l=0; l<BasisOnMeshType::dim(); l++)
        {
          const int delta_jl = (j == l? 1 : 0);
          const int Ii = delta_ik * delta_jl;
          
          const double Cc = inverseRightCauchyGreen[j][i] * rightCauchyGreen[l][k];
          const double Pp = (Ii - 1./BasisOnMeshType::dim() * Cc);
          
          pSbar += Pp * fictitiousPK2Stress[l][k];
          
          //cSbar += rightCauchyGreen[l][k] * fictitiousPK2Stress[l][k];
          //ccs += Cc * fictitiousPK2Stress[l][k];
        }
      }
      
      // isochoric stress
      const double sIso = factorJ23 * pSbar;
      pk2StressIsochoric[j][i] = sIso;
      
      // total stress is sum of volumetric and isochoric part
      pK2Stress[j][i] = sVol + sIso;
      
      //if (i == j)
      //  LOG(DEBUG) << "  ccs: " << ccs << " C:Sbar: " << cSbar << ", factorJ23: " << factorJ23 << ", Svol_" << i << j << " = " << sVol << ", Siso_" << i << j << " = " << sIso << ", S = " << pK2Stress[j][i]; 
    }
  }
  
  //check symmetry of PK2 stress
  const double errorTolerance = 1e-14;
  if (VLOG_IS_ON(2))
  {
    bool pK2IsSymmetric = true;
    for (int a=0; a<BasisOnMeshType::dim(); a++)
    {
      for (int b=0; b<BasisOnMeshType::dim(); b++)
      {
        if (fabs(pK2Stress[a][b] - pK2Stress[b][a]) > errorTolerance)
        {
          LOG(ERROR) << "pK2Stress["<<a<<"]["<<b<<"] != pK2Stress["<<b<<"]["<<a<<"] ("<<pK2Stress[b][a]<<" != "<<pK2Stress[a][b]<<")";
          pK2IsSymmetric = false;
        }
      }
    }
    if (pK2IsSymmetric)
      VLOG(2) << "PK2 stress tensor is symmetric!";
  }
  
  return pK2Stress;
}

template<typename BasisOnMeshType, typename Term>
void SolidMechanicsUtility<BasisOnMeshType, Term>:: 
checkFictitiousPK2Stress(const Tensor2<BasisOnMeshType::dim()> &fictitiousPK2Stress,
                         const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen, 
                         const double deformationGradientDeterminant,
                         const std::array<double,2> reducedInvariants)
{
  if (VLOG_IS_ON(1))
  {
    // explicit formula in Holzapfel p.249
   
    double factor23 = -2./3;
    if (BasisOnMeshType::dim() == 2)
    {
      factor23 = -1./2;
    }
    const double factorJ23 = pow(deformationGradientDeterminant, factor23);
    
    const double c0 = SEMT::Parameter<0>::get_value();
    const double c1 = SEMT::Parameter<1>::get_value();
    //const double c0 = PARAM(0).get_value();   //< material parameter
    //const double c1 = PARAM(1).get_value();   //< material parameter

    const double Ibar1 = reducedInvariants[0];
    
    const double gamma1 = 2*(c0 + c1*Ibar1);
    const double gamma2 = -2*c1;
    
    const double errorTolerance = 1e-14;
    bool mismatch = false;
    for (int a = 0; a < BasisOnMeshType::dim(); a++)
    {
      for (int b = 0; b < BasisOnMeshType::dim(); b++)
      {
        const int delta_ab = (a == b? 1 : 0);
        double cBar = factorJ23 * rightCauchyGreen[b][a];
        double sBar = gamma1*delta_ab + gamma2*cBar;
        
        if (fabs(sBar - fictitiousPK2Stress[b][a]) > errorTolerance)
        {
          LOG(ERROR) << "mismatch in Sbar_" << a << b << ": derived: " << fictitiousPK2Stress << ", explicit formula: " << sBar;
          mismatch = true;
        }
      }
    }
    if (!mismatch)
      VLOG(2) << "Sbar is correct!";
    
  }
}

template<typename BasisOnMeshType, typename Term>
Tensor2<BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeGreenLagrangeStrain(const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen)
{
  Tensor2<BasisOnMeshType::dim()> greenLagrangeStrain;
  
  for (int columnIndex = 0; columnIndex < BasisOnMeshType::dim(); columnIndex++)
  {
    for (int rowIndex = 0; rowIndex < BasisOnMeshType::dim(); rowIndex++)
    {
      const int delta_ij = (columnIndex == rowIndex? 1 : 0);
      greenLagrangeStrain[columnIndex][rowIndex] = 0.5*(rightCauchyGreen[columnIndex][rowIndex] - delta_ij);
    }
  }
  return greenLagrangeStrain;
}
  
template<typename BasisOnMeshType, typename Term>
ElasticityTensor SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeElasticityTensorCoupledStrainEnergy(const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen, 
                                           const Tensor2<BasisOnMeshType::dim()> &inverseRightCauchyGreen, 
                                           const std::array<double,3> invariants)
{
  // the 21 distinct indices (i,j,k,l) of different values of C_{ijkl}
  int indices[21][4] = {
    {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
    {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
    {2,2,2,2}
  };
  
  auto dPsidI1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunction, Term::I1);
  auto dPsidI2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunction, Term::I2);
  auto dPsidI3Expression = SEMT::deriv_t(Term::strainEnergyDensityFunction, Term::I3);
  
  auto d2PsidI1I1Expression = SEMT::deriv_t(dPsidI1Expression, Term::I1);
  auto d2PsidI1I2Expression = SEMT::deriv_t(dPsidI1Expression, Term::I2);
  auto d2PsidI1I3Expression = SEMT::deriv_t(dPsidI1Expression, Term::I3);
  auto d2PsidI2I2Expression = SEMT::deriv_t(dPsidI2Expression, Term::I2);
  auto d2PsidI2I3Expression = SEMT::deriv_t(dPsidI2Expression, Term::I3);
  auto d2PsidI3I3Expression = SEMT::deriv_t(dPsidI3Expression, Term::I3);
  
  std::vector<double> invariantsVector(invariants.begin(), invariants.end());
  
  //const double dPsidI1 = dPsidI1Expression.apply(invariantsVector);
  const double dPsidI2 = dPsidI2Expression.apply(invariantsVector);
  const double dPsidI3 = dPsidI3Expression.apply(invariantsVector);
  const double d2PsidI1I1 = d2PsidI1I1Expression.apply(invariantsVector);
  const double d2PsidI1I2 = d2PsidI1I2Expression.apply(invariantsVector);
  const double d2PsidI1I3 = d2PsidI1I3Expression.apply(invariantsVector);
  const double d2PsidI2I2 = d2PsidI2I2Expression.apply(invariantsVector);
  const double d2PsidI2I3 = d2PsidI2I3Expression.apply(invariantsVector);
  const double d2PsidI3I3 = d2PsidI3I3Expression.apply(invariantsVector);
  
  const double I1 = invariants[0];
  const double I3 = invariants[2];
  
  // formula cf. Holzapfel "Nonlinear Solid Mechanics" p.261
  // compute factors 
  const double factor1 = 4*(d2PsidI1I1 + 2*I1*d2PsidI1I2 + dPsidI2 + MathUtility::sqr(I1)*d2PsidI2I2);
  const double factor2 = -4*(d2PsidI1I2 + I1*d2PsidI2I2);
  const double factor3 = 4*(I3*d2PsidI1I3 + I1*I3*d2PsidI2I3);
  const double factor4 = 4*d2PsidI2I2;
  const double factor5 = -4*I3*d2PsidI2I3;
  const double factor6 = 4*(I3*dPsidI3 + MathUtility::sqr(I3)*d2PsidI3I3);
  const double factor7 = -4*I3*dPsidI3;
  const double factor8 = -4*dPsidI2;
  
  std::array<double,21> elasticity({0});
  // loop over distinct entries in elasticity tensor
  for (int entryNo = 0; entryNo < 21; entryNo++)
  {
    // rename indices of current entry
    const int i = indices[entryNo][0];
    const int j = indices[entryNo][1];
    const int k = indices[entryNo][2];
    const int l = indices[entryNo][3];
    
    // note that rightCauchyGreen and inverseRightCauchyGreen are stored column-major i.e. C_{ij} = rightCauchyGreen[j][i]
    
    // factor1 * (I dyad I)
    elasticity[entryNo] += factor1 * (i==j) * (k==l);   
    
    // factor2 * (I dyad C + C dyad I)
    elasticity[entryNo] += factor2 * ((i==j) * rightCauchyGreen[l][k] + rightCauchyGreen[j][i] * (k==l));  
    
    // factor3 * (I dyad C^-1 + C^-1 dyad I)
    elasticity[entryNo] += factor3 * ((i==j) * inverseRightCauchyGreen[l][k] + inverseRightCauchyGreen[j][i] * (k==l));  
    
    // factor4 * (C dyad C)
    elasticity[entryNo] += factor4 * (rightCauchyGreen[j][i] * rightCauchyGreen[l][k]);
    
    // factor5 * (C dyad C^-1 + C^-1 dyad C)
    elasticity[entryNo] += factor5 * (rightCauchyGreen[j][i] * inverseRightCauchyGreen[l][k] + inverseRightCauchyGreen[j][i] * rightCauchyGreen[l][k]);
    
    // factor6 * (C^-1 dyad C^-1)
    elasticity[entryNo] += factor6 * (inverseRightCauchyGreen[j][i] * inverseRightCauchyGreen[l][k]);
    
    // factor7 * (C^-1 odot C^-1), where A odot A := 1/2*(A_{ik}*A_{jl} + A_{il}*A_{jk})
    elasticity[entryNo] += factor7 
      * 1./2*(inverseRightCauchyGreen[k][i]*inverseRightCauchyGreen[l][j] + inverseRightCauchyGreen[l][i]*inverseRightCauchyGreen[k][j]);
      
    // factor8 * II (where II = delta_{ik}*delta_{jl}*e_i dyad e_j dyad e_k dyad e_l
    elasticity[entryNo] += factor8 * (i==k) * (j==l);
  }
  
  return elasticity;
}

template<typename BasisOnMeshType, typename Term>
double SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeElasticityTensorEntry(const int i, const int j, const int k, const int l,const double pressure, 
                        const double pressureTilde,
                        const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen,
                        const Tensor2<BasisOnMeshType::dim()> &inverseRightCauchyGreen,
                        const Tensor2<BasisOnMeshType::dim()> &fictitiousPK2Stress,
                        const Tensor2<BasisOnMeshType::dim()> &pk2StressIsochoric,
                        const double deformationGradientDeterminant,
                        const std::array<double,2> reducedInvariants)
{
  auto dPsi_dIbar1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar1);
  auto dPsi_dIbar2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunctionIsochoric, Term::Ibar2);
  
  auto d2Psi_dIbar1Ibar1Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar1);
  auto d2Psi_dIbar1Ibar2Expression = SEMT::deriv_t(dPsi_dIbar1Expression, Term::Ibar2);
  auto d2Psi_dIbar2Ibar2Expression = SEMT::deriv_t(dPsi_dIbar2Expression, Term::Ibar2);
  
  std::vector<double> reducedInvariantsVector(reducedInvariants.begin(), reducedInvariants.end());
  
  //const double dPsi_dIbar1 = dPsi_dIbar1Expression.apply(reducedInvariantsVector);
  const double dPsi_dIbar2 = dPsi_dIbar2Expression.apply(reducedInvariantsVector);
  const double d2Psi_dIbar1Ibar1 = d2Psi_dIbar1Ibar1Expression.apply(reducedInvariantsVector);
  const double d2Psi_dIbar1Ibar2 = d2Psi_dIbar1Ibar2Expression.apply(reducedInvariantsVector);
  const double d2Psi_dIbar2Ibar2 = d2Psi_dIbar2Ibar2Expression.apply(reducedInvariantsVector);
  
  const double Ibar1 = reducedInvariants[0];
  
  // formula for C_iso: Holzapfel "Nonlinear Solid Mechanics" p.255
  // formula for C_bar: Holzapfel "Nonlinear Solid Mechanics" p.262
  // compute factors for Cbar
  const double factor1 = 4*(d2Psi_dIbar1Ibar1 + 2*Ibar1*d2Psi_dIbar1Ibar2 + dPsi_dIbar2 + MathUtility::sqr(Ibar1)*d2Psi_dIbar2Ibar2);
  const double factor2 = -4*(d2Psi_dIbar1Ibar2 + Ibar1*d2Psi_dIbar2Ibar2);
  const double factor3 = 4*d2Psi_dIbar2Ibar2;
  const double factor4 = -4*dPsi_dIbar2;
  
  const double J = deformationGradientDeterminant;
  
  double factor23 = -2./3;
  double factor43 = -4./3;
  if (BasisOnMeshType::dim() == 2)
  {
    factor23 = -1./2;
    factor43 = -1;
  }
  const double factorJ23 = pow(J,factor23);
  const double factorJ43 = pow(J,factor43);
  
  // Cbar = J^{-2/3}*C
  
  // P = II - 1/3 C^-1 dyad C    (4th order tensor, p.229)
  // II_ijkl = delta_ik*delta_jl (p.23)
  // (P^T)_ijkl = P_klij
  
#if 0
  //check symmetry of right Cauchy-Green tensor
  bool cIsSymmetric = true;
  for (int a=0; a<3; a++)
  {
    for (int b=0; b<3; b++)
    {
      if (rightCauchyGreen[a][b] != rightCauchyGreen[b][a])
      {
        std::cout << "rightCauchyGreen["<<a<<"]["<<b<<"] != " << "rightCauchyGreen["<<b<<"]["<<a<<"] ("<<rightCauchyGreen[b][a]<<" != "<<rightCauchyGreen[a][b]<<")"<<std::endl;
        cIsSymmetric = false;
      }
    }
  }
  if (cIsSymmetric)
    LOG(DEBUG) << "right cauchy green tensor is symmetric!";
  
  double Cbar[3][3][3][3];
  double Ii1[3][3][3][3];
  double Ii2[3][3][3][3];
  double Icci[3][3][3][3];
  double Cc[3][3][3][3];
  
  std::stringstream s;
  for (int a=0; a<3; a++)
  {
    for (int b=0; b<3; b++)
    {
      for (int c=0; c<3; c++)
      {
        for (int d=0; d<3; d++)
        {
          // C = J^(2/3) * Cbar  =>   Cbar = J^(-2/3) * C = factorJ23*rightCauchyGreen
          double J43Cbar = 0.;  // J43Cbar = J^{4/3}*Cbar_abcd => J^{-4/3} * J43Cbar = Cbar_abcd
          
          // factor1 * (I dyad I)
          J43Cbar += factor1 * (a==b) * (c==d);   
          
          // factor2 * (I dyad Cbar + Cbar dyad I)
          J43Cbar += factor2 * factorJ23 * ((a==b) * rightCauchyGreen[d][c] + rightCauchyGreen[b][a] * (c==d));  
          
          // factor3 * (Cbar dyad Cbar)
          J43Cbar += factor3 * (factorJ23 * rightCauchyGreen[b][a] * factorJ23 * rightCauchyGreen[d][c]);
          
          // factor4 * II (where II = delta_{ac}*delta_{bd}*e_a dyad e_b dyad e_c dyad e_d
          J43Cbar += factor4 * (a==c) * (b==d);
          
          const double Cbar_abcd = factorJ43 * J43Cbar;
          Cbar[a][b][c][d] = Cbar_abcd;
          
          Ii1[a][b][c][d] = (a==b? 1 : 0) * (c==d? 1 : 0);
          Ii2[a][b][c][d] = (a==c? 1 : 0) * (b==d? 1 : 0);
          Icci[a][b][c][d] = (a==b? 1 : 0) * factorJ23 * rightCauchyGreen[d][c] + factorJ23 * rightCauchyGreen[b][a] * (c==d? 1 : 0);
          Cc[a][b][c][d] = factorJ23 * rightCauchyGreen[b][a] * factorJ23 * rightCauchyGreen[d][c];
          
          s << Cbar_abcd << " ";
        }
      }
      s << std::endl;
    }
  }
  
  LOG(DEBUG) << "    Cbar_abcd: " << std::endl << s.str();
  checkSymmetry(Cbar,"Cbar");
  checkSymmetry(Ii1,"Ii1");
  checkSymmetry(Ii2,"Ii2");
  checkSymmetry(Icci,"Icci");
  checkSymmetry(Cc,"Cc");
  
  
  double CbarIiT[3][3][3][3];
  s.str("");
  for (int a=0; a<3; a++)
  {
    for (int b=0; b<3; b++)
    {
      for (int k=0; k<3; k++)
      {
        for (int l=0; l<3; l++)
        {
          double CbarIiT_abkl = 0;
          
          for (int c=0; c<3; c++)
          {
            for (int d=0; d<3; d++)
            {
              double J43Cbar = 0.;  // J43Cbar = J^{4/3}*Cbar_abcd => J^{-4/3} * J43Cbar = Cbar_abcd
              
              // factor1 * (I dyad I)
              J43Cbar += factor1 * (a==b) * (c==d);   
              
              // factor2 * (I dyad Cbar + Cbar dyad I)
              J43Cbar += factor2 * factorJ23 * ((a==b) * rightCauchyGreen[d][c] + rightCauchyGreen[b][a] * (c==d));  
              
              // factor3 * (Cbar dyad Cbar)
              J43Cbar += factor3 * (factorJ23 * rightCauchyGreen[b][a] * factorJ23 * rightCauchyGreen[d][c]);
              
              // factor4 * II (where II = delta_{ac}*delta_{bd}*e_a dyad e_b dyad e_c dyad e_d
              J43Cbar += factor4 * (a==c) * (b==d);
              
              const double Cbar_abcd = factorJ43 * J43Cbar;
              
              const int delta_kc = (k == c? 1 : 0);
              const int delta_ld = (l == d? 1 : 0);
          
        
              const int Ii_klcd = delta_kc * delta_ld;
              const int IiT_cdkl = Ii_klcd;
              
              CbarIiT_abkl += Cbar_abcd * IiT_cdkl;
            }
          }
        
          CbarIiT[a][b][k][l] = CbarIiT_abkl;
          s << CbarIiT_abkl << " ";
        }
      }
      s << std::endl;
    }
  }
  
  LOG(DEBUG) << "    CbarIiT_abkl: " << std::endl << s.str();
  checkSymmetry(CbarIiT,"CbarIiT");
#endif  
  
  // ----------- Ciso (Holzapfel p.255) -------------------------
  //                          ij ab    cd    kl
  // compute contribution from P : Cbar : P^T
  double PCbarPT_ijkl = 0.;
  
  // row index
  for (int a=0; a<BasisOnMeshType::dim(); a++)
  {
    // column index
    for (int b=0; b<BasisOnMeshType::dim(); b++)
    {
      const int delta_ia = (i == a? 1 : 0);
      const int delta_jb = (j == b? 1 : 0);
      
      // P = II - 1/3 C^-1 dyad C    (4th order tensor, p.229)
      const int Ii_ijab = delta_ia * delta_jb;
      const double P_ijab = Ii_ijab - 1./BasisOnMeshType::dim() * inverseRightCauchyGreen[j][i] * rightCauchyGreen[b][a];
      //const double P_ijab = Ii_ijab;
     
      double CbarPT_abkl = 0.0;
      
      // row index
      for (int c=0; c<BasisOnMeshType::dim(); c++)
      {
        const int delta_kc = (k == c? 1 : 0);
        
        // column index
        for (int d=0; d<BasisOnMeshType::dim(); d++)
        {
          const int delta_ld = (l == d? 1 : 0);
          
          const int Ii_klcd = delta_kc * delta_ld;
          const double P_klcd = Ii_klcd - 1./BasisOnMeshType::dim() * inverseRightCauchyGreen[l][k] * rightCauchyGreen[d][c];
          //const double P_klcd = Ii_klcd;
          
          const double PT_cdkl = P_klcd;  // Holzapfel p.23
           
          double J43Ccbar_abcd = 0.;  // J43Cbar = J^{4/3}*Cbar_abcd => J^{-4/3} * J43Cbar = Cbar_abcd
          // Cbar: Holzapfel p.262
          
          const int delta_ab = (a == b? 1 : 0);
          const int delta_cd = (c == d? 1 : 0);
          const int delta_ac = (a == c? 1 : 0);
          const int delta_bd = (b == d? 1 : 0);
          
          // factor1 * (I dyad I)
          J43Ccbar_abcd += factor1 * delta_ab * delta_cd;   
          
          // factor2 * (I dyad Cbar + Cbar dyad I)
          J43Ccbar_abcd += factor2 * factorJ23 * (delta_ab * rightCauchyGreen[d][c] + rightCauchyGreen[b][a] * delta_cd);  
          
          // factor3 * (Cbar dyad Cbar)
          J43Ccbar_abcd += factor3 * (factorJ23 * rightCauchyGreen[b][a] * factorJ23 * rightCauchyGreen[d][c]);
          
          // factor4 * II (where II = delta_{ac}*delta_{bd}*e_a dyad e_b dyad e_c dyad e_d
          J43Ccbar_abcd += factor4 * delta_ac * delta_bd;
          
          const double Cbar_abcd = factorJ43 * J43Ccbar_abcd;
          
          //LOG(DEBUG) << "    Cbar_abcd [" << a << b << c << d << "]: " << Cbar_abcd << ", add " << Cbar_abcd << "*" << PT_cdkl << "=" << Cbar_abcd * PT_cdkl;
          
          CbarPT_abkl += Cbar_abcd * PT_cdkl;
        }
      }
      
      PCbarPT_ijkl += P_ijab * CbarPT_abkl;
      
      //LOG(DEBUG) << "   CbarPT_abkl [" << a << b << k << l <<"]: " << CbarPT_abkl;
      
    }
  }

  //LOG(DEBUG) << "   PCbarPT_ijkl [" << i << j << k << l << "]: " << PCbarPT_ijkl;
  
  // compute Tr(Sbar)
  // Sbar = fictitiousPK2Stress, Tr(•) = (•):C
  double TrSbar = 0.0;
  for (int a=0; a<BasisOnMeshType::dim(); a++)
  {
    // column index
    for (int b=0; b<BasisOnMeshType::dim(); b++)
    {
      TrSbar += fictitiousPK2Stress[b][a] * rightCauchyGreen[b][a];
    }
  }
  
  // Holzapfel p.255
  const double CcTerm1_ijkl = 1./2*(inverseRightCauchyGreen[k][i]*inverseRightCauchyGreen[l][j] + inverseRightCauchyGreen[l][i]*inverseRightCauchyGreen[k][j]);
  const double CcTerm2_ijkl = inverseRightCauchyGreen[j][i]*inverseRightCauchyGreen[l][k];
  const double pTilde_ijkl = CcTerm1_ijkl - 1./BasisOnMeshType::dim() * CcTerm2_ijkl;
  
  const double trTerm_ijkl = -factor23 * factorJ23 * TrSbar * pTilde_ijkl;
  
  // compute last term -2/3(C^{-1} dyad S_iso + S_iso dyad C^{-1})
  const double lastTerm_ijkl = factor23*(inverseRightCauchyGreen[j][i]*pk2StressIsochoric[l][k] + pk2StressIsochoric[j][i]*inverseRightCauchyGreen[l][k]);
  
  const double Ciso = PCbarPT_ijkl + trTerm_ijkl + lastTerm_ijkl;
  
  // ----------- Cvol p.254 -------------------------
  const double Cvol = J*pressureTilde*CcTerm2_ijkl - 2*J*pressure*CcTerm1_ijkl;
  
  //LOG(DEBUG) << "    Cvol: " << Cvol << ", Ciso: " << Ciso << " = " << PCbarPT_ijkl << " + " << trTerm_ijkl << " + " << lastTerm_ijkl;
  
  return Cvol + Ciso;
}

template<typename BasisOnMeshType, typename Term>
ElasticityTensor SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computeElasticityTensor(const double pressure, 
                        const double pressureTilde,
                        const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen,
                        const Tensor2<BasisOnMeshType::dim()> &inverseRightCauchyGreen,
                        const Tensor2<BasisOnMeshType::dim()> &fictitiousPK2Stress,
                        const Tensor2<BasisOnMeshType::dim()> &pk2StressIsochoric,
                        const double deformationGradientDeterminant,
                        const std::array<double,2> reducedInvariants)
{
  // the 21 distinct indices (i,j,k,l) of different values of C_{ijkl}
  // the symmetries are: C_ijkl = C_klij, C_ijkl = C_jikl = C_ijlk = C_jilk
  int indices[21][4] = {
    {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
    {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
    {2,2,2,2}
  };
  
  ElasticityTensor elasticity;
  // loop over distinct entries in elasticity tensor
  for (int entryNo = 0; entryNo < 21; entryNo++)
  {
    // rename indices of current entry
    const int i = indices[entryNo][0];
    const int j = indices[entryNo][1];
    const int k = indices[entryNo][2];
    const int l = indices[entryNo][3];
    
    VLOG(2) << " ----";
    VLOG(2) << " elasticity ijkl = " << i << j << k << l << " (entry " << entryNo << ")";
    double entry_ijkl = computeElasticityTensorEntry(i,j,k,l,pressure,pressureTilde,rightCauchyGreen,inverseRightCauchyGreen,fictitiousPK2Stress,pk2StressIsochoric,deformationGradientDeterminant,reducedInvariants);
    
    //const double tol = 1e-15;
    VLOG(2) << " ijkl: " << entry_ijkl;
    
    //double entry_jikl = computeElasticityTensorEntry(j,i,k,l,pressure,pressureTilde,rightCauchyGreen,inverseRightCauchyGreen,fictitiousPK2Stress,pk2StressIsochoric,deformationGradientDeterminant,reducedInvariants);
    //VLOG(2) << " jikl: " << entry_jikl << " (diff: " << entry_jikl-entry_ijkl << ")";
    
    //double entry_ijlk = computeElasticityTensorEntry(i,j,l,k,pressure,pressureTilde,rightCauchyGreen,inverseRightCauchyGreen,fictitiousPK2Stress,pk2StressIsochoric,deformationGradientDeterminant,reducedInvariants);
    //VLOG(2) << " ijlk: " << entry_ijlk << " (diff: " << entry_ijlk-entry_ijkl << ")";
    
    //double entry_klij = computeElasticityTensorEntry(k,l,i,j,pressure,pressureTilde,rightCauchyGreen,inverseRightCauchyGreen,fictitiousPK2Stress,pk2StressIsochoric,deformationGradientDeterminant,reducedInvariants);
    //VLOG(2) << " klij: " << entry_klij << " (diff: " << entry_klij-entry_ijkl << ")";
    
    //assert(fabs(entry_ijkl - entry_jikl) < tol);
    //assert(fabs(entry_ijkl - entry_ijlk) < tol);
    //assert(fabs(entry_ijkl - entry_klij) < tol);
    
    elasticity[entryNo] = entry_ijkl;
  }
  
  return elasticity;
}

template<typename BasisOnMeshType, typename Term>
double SolidMechanicsUtility<BasisOnMeshType, Term>:: 
computePressureFromDisplacements(double deformationGradientDeterminant, 
                                 const Tensor2<BasisOnMeshType::dim()> &rightCauchyGreen, 
                                 const Tensor2<BasisOnMeshType::dim()> &PK2Stress)
{
  double pressure = 0;
  double factor = 0;
  
  // row index
  for (int i=0; i<BasisOnMeshType::dim(); i++)
  {
    // column index
    for (int j=0; j<BasisOnMeshType::dim(); j++)
    {
      factor += rightCauchyGreen[j][i] * PK2Stress[j][i];
    }
  }
  pressure = -1.0 / (BasisOnMeshType::dim()*deformationGradientDeterminant) * factor;
  
  return pressure;
}
 
};  // namespace