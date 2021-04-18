#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available

#include "equation/mooney_rivlin_incompressible.h"
#include "utility/math_utility.h"

namespace SpatialDiscretization
{

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
std::array<double_v_t,5> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
std::array<double_v_t,5> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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
#ifndef NDEBUG
  else
  {
    reducedInvariants[2] = 1;
    reducedInvariants[3] = 0;
    reducedInvariants[4] = 0;
  }
#endif

  if (Vc::any_of(deformationGradientDeterminant <= 0))
  {
#ifdef HAVE_STDSIMD
    LOG(ERROR) << "J=det F is negative. Result will be unphysical.\n"
      << "For dynamic problems, reduce time step width, for static problems, add smaller \"loadFactors\" or reduce load.";
#else
    LOG(ERROR) << "J=det F is negative: " << deformationGradientDeterminant << ". Result will be unphysical.\n"
      << "For dynamic problems, reduce time step width, for static problems, add smaller \"loadFactors\" or reduce load.";
#endif

    Vc::where(deformationGradientDeterminant <= 0, reducedInvariants[0]) = 3;
    Vc::where(deformationGradientDeterminant <= 0, reducedInvariants[1]) = 0;
    Vc::where(deformationGradientDeterminant <= 0, reducedInvariants[3]) = 3;
    Vc::where(deformationGradientDeterminant <= 0, reducedInvariants[4]) = 0;
  }

  return reducedInvariants;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
template<typename double_v_t>
Tensor2<3,double_v_t> HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
computePSbar(const Tensor2<3,double_v_t> &fictitiousPK2Stress, const Tensor2<3,double_v_t> &rightCauchyGreen)
{
  // only needed for debugging in materialTesting
  double_v_t determinant;
  double_v_t approximateMeshWidth{0};
  Tensor2<3,double_v_t> inverseRightCauchyGreen = MathUtility::computeInverse<3>(rightCauchyGreen, approximateMeshWidth, determinant);

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

