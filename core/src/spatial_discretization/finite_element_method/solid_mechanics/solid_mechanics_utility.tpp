#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

namespace SpatialDiscretization
{  

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<Vec3,BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Term>:: 
computeDeformationGradient(const std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &geometryReference,
                           const std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &displacement, 
                           const std::array<Vec3,BasisOnMeshType::dim()> &jacobian, 
                           const std::array<double, BasisOnMeshType::dim()> xi)
{
 
  // compute the deformation gradient x_i,j = d_ij + u_i,j
  // where j is dimensionColumn and i i component of the used Vec3's
  std::array<Vec3,BasisOnMeshType::dim()> deformationGradient;
 
  const int nDofsPerElement = BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement();
  
  // loop over dimension, i.e. columns of deformation gradient
  for (int dimensionColumn = 0; dimensionColumn < BasisOnMeshType::dim(); dimensionColumn++)
  {
    Vec3 dudxi({0});   // handle full-dimension vector of displacement (i.e. (x,y,z))
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      double dphi_dxi = BasisOnMeshType::HighOrderBasisOnMesh::dphi_dxi(dofIndex, dimensionColumn, xi);
      dudxi += dphi_dxi * displacement[dofIndex];   // vector-valued addition
    }
    
    // multiply du/dxi with dxi/dx to obtain du/dx
    deformationGradient[dimensionColumn] = dudxi * jacobian[dimensionColumn];
    
    // add Kronecker delta to obtain x_i,j = delta_ij + u_i,j
    deformationGradient[dimensionColumn][dimensionColumn] += 1;
  }
}

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<Vec3,BasisOnMeshType::dim()> SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Term>:: 
computeRightCauchyGreenTensor(const std::array<Vec3,BasisOnMeshType::dim()> &deformationGradient)
{
  // compute C = F^T*F where F is the deformationGradient and C is the right Cauchy-Green Tensor
  // the quantities are 3x3 tensors for the 3D case (3x2 for 2D, but not sure if that is required, 2D is not used anyway)
  std::array<Vec3,BasisOnMeshType::dim()> rightCauchyGreenTensor({Vec3{0}});
  
  // loop over dimension, i.e. columns of right cauchy green tensor
  for (int dimensionColumn = 0; dimensionColumn < BasisOnMeshType::dim(); dimensionColumn++)
  {
    // loop over row of tensor
    for (int dimensionRow = 0; dimensionRow < 3; dimensionRow++)
    {
      for (int k = 0; k < 3; k++)
      {
        rightCauchyGreenTensor[dimensionColumn][dimensionRow] += 
          deformationGradient[dimensionColumn][k] * deformationGradient[k][dimensionRow];
      }
    }
  }
  
  return rightCauchyGreenTensor;
}

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<double,3> SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Term>:: 
computeInvariants(const std::array<Vec3,BasisOnMeshType::dim()> &rightCauchyGreen, const double determinant)
{
  std::array<double,3> invariants;
  
  // I1 = tr(C)
  invariants[0] = rightCauchyGreen[0][0] + rightCauchyGreen[1][1] + rightCauchyGreen[2][2];
  
  // tr(C^2)
  double traceC2 = 0;
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
  invariants[2] = determinant;
}

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<Vec3,3> SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Term>:: 
computePK2Stress(const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &inverseRightCauchyGreen, const std::array<double,3> invariants)
{
  // compute the PK2 stress tensor as S=2*dPsi/dC
  // for explanation see pdf document
  auto dPsidI1Expression = SEMT::deriv_t(Term::strainEnergyDensityFunction, Term::I1);
  auto dPsidI2Expression = SEMT::deriv_t(Term::strainEnergyDensityFunction, Term::I2);
  auto dPsidI3Expression = SEMT::deriv_t(Term::strainEnergyDensityFunction, Term::I3);
  
  std::vector<double> invariantsVector(invariants.begin(), invariants.end());
  
  const double dPsidI1 = dPsidI1Expression.apply(invariantsVector);
  const double dPsidI2 = dPsidI2Expression.apply(invariantsVector);
  const double dPsidI3 = dPsidI3Expression.apply(invariantsVector);
  
  const double I1 = invariants[0];
  const double I3 = invariants[2];
  
  double factor1 = dPsidI1 + I1*dPsidI2;
  double factor2 = dPsidI2;
  double factor3 = I3*dPsidI3;
  
  std::array<Vec3,3> pK2Stress;
  
  // row index
  for (int i=0; i<3; i++)
  {
    // column index
    for (int j=0; j<3; j++)
    {
      pK2Stress[j][i] = 2*(factor1 * (i == j? 1 : 0) - factor2 * rightCauchyGreen[j][i] + factor3 * inverseRightCauchyGreen[j][i]);
    }
  }
  
  return pK2Stress;
}

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<double,21> SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Term>:: 
computeElasticityTensor(const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &inverseRightCauchyGreen, const std::array<double,3> invariants)
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

template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
double SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Term>:: 
computePressureFromDisplacements(double deformationGradientDeterminant, const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &PK2Stress)
{
  double pressure = 0;
  double factor = 0;
  
  // row index
  for (int i=0; i<3; i++)
  {
    // column index
    for (int j=0; j<3; j++)
    {
      factor += rightCauchyGreen[j][i] * PK2Stress[j][i];
    }
  }
  pressure = -1.0 / (3*deformationGradientDeterminant) * factor;
}
 
};  // namespace