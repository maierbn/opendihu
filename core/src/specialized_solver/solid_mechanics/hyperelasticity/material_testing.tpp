#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "equation/mooney_rivlin_incompressible.h"

namespace SpatialDiscretization
{

template<typename Term,typename MeshType, int nDisplacementComponents>
double HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
computeSbarC(const Tensor2<3> &Sbar, const Tensor2<3> &C)
{
  double SbarC = 0;
  for (int a = 0; a < 3; a++)
  {
    for (int b = 0; b < 3; b++)
    {
      SbarC += Sbar[b][a] * C[b][a];
    }
  }
  return SbarC;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
materialTesting(const double pressure,                           //< [in] pressure value p
                const Tensor2<3> &rightCauchyGreen,                //< [in] C
                const Tensor2<3> &inverseRightCauchyGreen,         //< [in] C^{-1}
                const std::array<double,5> reducedInvariants,      //< [in] the reduced invariants Ibar_1, Ibar_2
                const double deformationGradientDeterminant,       //< [in] J = det(F)
                Vec3 fiberDirection,                               //< [in] a0, direction of fibers
                Tensor2<3> &fictitiousPK2Stress,                   //< [in] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                Tensor2<3> &pk2StressIsochoric                    //< [in] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
               )
{
#if 0
  LOG(DEBUG) << "materialTesting, parameters: " << this->materialParameters_ << ", C: " << rightCauchyGreen;

  const int D = 3;

  // set reference values
  double JReference                  = deformationGradientDeterminant;
  Tensor2<3> CReference              = rightCauchyGreen;
  Tensor2<3> CInvReference           = inverseRightCauchyGreen;
  std::array<double,5> IbarReference = reducedInvariants;

  // do not do testing for reference configuration
  if (fabs(JReference-1) < 1e-4)
    return;

  // compute CBarReference = JReference^{-2/3} * CReference
  Tensor2<3> CBarReference;
  for (int a = 0; a < 3; a++)
  {
    for (int b = 0; b < 3; b++)
    {
      CBarReference[b][a] = pow(JReference, -2./3) * CReference[b][a];
    }
  }
  double determinant;
  Tensor2<3> CBarInvReference = MathUtility::computeInverse<3>(CBarReference, determinant);


  // compute Sbar from input values
  Tensor2<3> SBarReference;
  Tensor2<3> SIsoReference;
  this->computePK2Stress(pressure, CReference, CInvReference, IbarReference, JReference, fiberDirection,
                                                  SBarReference, SIsoReference
                                                );

#if 0

  // compute Psi
  auto PsiExpression = Term::strainEnergyDensityFunctionIsochoric;
  std::vector<double> reducedInvariantsVector(IbarReference.begin(), IbarReference.end());
  const double PsiReference = PsiExpression.apply(reducedInvariantsVector);


  LOG(DEBUG) << "reference values: ";
  LOG(DEBUG) << "  J: " << JReference << ", C: " << CReference << ", CInv: " << CInvReference;
  LOG(DEBUG) << "  Ibar: " << IbarReference[0] << "," << IbarReference[1] << ", CBar: " << CBarReference;
  LOG(DEBUG) << "  Psi: " << PsiReference;
  LOG(DEBUG) << "  Sbar: " << SBarReference;
  LOG(DEBUG) << "  Siso: " << SIsoReference;

  // test if Sbar = 2 ∂ψ/∂Cbar by evaluating the derivative numerically (finite differences)
  Tensor2<3> SPsi;

  double h = 1e-8;

  if (true)
  {
    // loop over entries of Cbar
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb Cbar
        Tensor2<3> Cbar = CBarReference;
        Cbar[b][a] += h/2.;
        Cbar[a][b] += h/2.;

        // compute Ibar
        std::array<double,5> Ibar;
        Ibar[0] = Cbar[0][0] + Cbar[1][1] + Cbar[2][2];

        double trC2 = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            trC2 += Cbar[i][j] * Cbar[j][i];
          }
        }
        Ibar[1] = 0.5*((Ibar[0]*Ibar[0]) - trC2);
        Ibar[2] = 1;

        std::vector<double> reducedInvariantsVector(Ibar.begin(), Ibar.end());
        double Psi = PsiExpression.apply(reducedInvariantsVector);

        SPsi[b][a] = 2*(Psi - PsiReference) / h;

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double derivative = (Cbar[d][c] - CBarReference[d][c]) / h;

            int delta_ac = (a == c? 1: 0);
            int delta_bd = (b == d? 1: 0);
            int delta_ad = (a == d? 1: 0);
            int delta_bc = (b == c? 1: 0);

            // abcd
            double ii = delta_ac * delta_bd;
            double ss = 0.5*(ii + delta_ad * delta_bc);

            LOG(DEBUG) << a << b << c << d << " derivatives: num: " << derivative << ", ii: " << ii << ", ss: " << ss;
          }
        }
        //LOG(DEBUG) << "  Cbar: " << Cbar << ", Ibar: " << Ibar[0] << "," << Ibar[1] << ", Psi: " << Psi;
      }  // b
    }  // a

    // compare SPsi and SBarReference
    LOG(DEBUG) << "SPsi:          " << SPsi;
    LOG(DEBUG) << "SBarReference: " << SBarReference;

    double error = 0;

    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        error += MathUtility::sqr(SPsi[b][a] - SBarReference[b][a]);
      }
    }
    LOG(DEBUG) << "Sbar error: " << error/9.;
    LOG(DEBUG) << "--";

    //if (error > 1)
    //  LOG(ERROR) << "Sbar error is too high";
  }

  // -----------------------------------

  if (true)
  {
    // compute reference C
    Tensor4<D> CC;
    Tensor4<D> CCbarReference;
    Tensor4<3> CCIsoReference;
    computeElasticityTensor(CReference, inverseRightCauchyGreen, JReference, pressure, IbarReference, SBarReference, SIsoReference, fiberDirection,
                            CCbarReference, CCIsoReference, CC);

    Tensor4<D> CCbar;

    // test if CCbar = J^{-4/3} * 2 * ∂Sbar/∂Cbar by evaluating the 1st derivative numerically
    // loop over entries of Cbar
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb Cbar
        Tensor2<3> Cbar = CBarReference;
        Cbar[b][a] += h/2.;
        Cbar[a][b] += h/2.;

        // compute Ibar
        std::array<double,5> Ibar;
        Ibar[0] = Cbar[0][0] + Cbar[1][1] + Cbar[2][2];

        double trC2 = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            trC2 += Cbar[i][j] * Cbar[j][i];
          }
        }
        Ibar[1] = 0.5*((Ibar[0]*Ibar[0]) - trC2);
        Ibar[2] = 1;

        // I4 = a0 • C a0;
        double a0Ca0 = 0;
        for (int i = 0; i < 3; i++)
        {
          double ca0_i = 0;
          for (int j = 0; j < 3; j++)
          {
            ca0_i += Cbar[j][i] * fiberDirection[j];
          }
          a0Ca0 += fiberDirection[i] * ca0_i;
        }

        Ibar[3] = a0Ca0;

        // compute Cbar^-1
        Tensor2<3> CbarInv = MathUtility::computeInverse<3>(Cbar, determinant);

        // compute Sbar, using J=1 such that input C = Cbar
        Tensor2<3> SBar;
        Tensor2<3> Siso;
        this->computePK2Stress(pressure, Cbar, CbarInv, Ibar, 1, fiberDirection,
                                                        SBar, Siso);

        //LOG(DEBUG) << "  Cbar: " << Cbar << ", Ibar: " << Ibar[0] << "," << Ibar[1] << ", SBar: " << SBar;

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {

            double derivative = (SBar[d][c] - SBarReference[d][c]) / h * 2 * pow(JReference, -4/3.);
            CCbar[b][a][d][c] = derivative;

            LOG(DEBUG) << a << b << c << d << " SBar[d][c]: " << SBar[d][c] << " (diff:" << SBar[d][c]-SBarReference[d][c]
              << ") -> CCbar from S: " << derivative << "(ref: " << CCbarReference[d][c][b][a] << ")";
          }  // d
        }  // c
      }  // b
    }  // a


    // test if CCbar = J^{-4/3} * 4 * ∂^2Sbar/∂Cbar^2 by evaluating the 2nd derivative numerically

    Tensor4<D> CCbar2;
    if (false)
    {
      //h = 1e-10;

      // loop over entries of Cbar
      for (int a = 0; a < 3; a++)
      {
        for (int b = 0; b < 3; b++)
        {
          for (int c = 0; c < 3; c++)
          {
            for (int d = 0; d < 3; d++)
            {
              // evaluation CC(C_ab+h,C_cd+h)
              // perturb Cbar
              Tensor2<3> Cbar = CBarReference;
              Cbar[b][a] += h;
              Cbar[a][b] += h;

              Cbar[d][c] += h;
              Cbar[c][d] += h;

              // compute Ibar
              std::array<double,5> Ibar;
              Ibar[0] = Cbar[0][0] + Cbar[1][1] + Cbar[2][2];

              double trC2 = 0;
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  trC2 += Cbar[i][j] * Cbar[j][i];
                }
              }
              Ibar[1] = 0.5*((Ibar[0]*Ibar[0]) - trC2);
              Ibar[2] = 1;

              // I4 = a0 • C a0;
              double a0Ca0 = 0;
              for (int i = 0; i < 3; i++)
              {
                double ca0_i = 0;
                for (int j = 0; j < 3; j++)
                {
                  ca0_i += Cbar[j][i] * fiberDirection[j];
                }
                a0Ca0 += fiberDirection[i] * ca0_i;
              }

              Ibar[3] = a0Ca0;

              std::vector<double> reducedInvariantsVector(Ibar.begin(), Ibar.end());
              double Psi11 = PsiExpression.apply(reducedInvariantsVector);

              // evaluation CC(C_ab+h,C_cd)
              // perturb Cbar
              Cbar = CBarReference;
              Cbar[b][a] += h;
              Cbar[a][b] += h;

              // compute Ibar
              Ibar[0] = Cbar[0][0] + Cbar[1][1] + Cbar[2][2];

              trC2 = 0;
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  trC2 += Cbar[i][j] * Cbar[j][i];
                }
              }
              Ibar[1] = 0.5*((Ibar[0]*Ibar[0]) - trC2);
              Ibar[2] = 1;

              // I4 = a0 • C a0;
              a0Ca0 = 0;
              for (int i = 0; i < 3; i++)
              {
                double ca0_i = 0;
                for (int j = 0; j < 3; j++)
                {
                  ca0_i += Cbar[j][i] * fiberDirection[j];
                }
                a0Ca0 += fiberDirection[i] * ca0_i;
              }

              Ibar[3] = a0Ca0;

              std::vector<double> reducedInvariantsVector2(Ibar.begin(), Ibar.end());
              double Psi10 = PsiExpression.apply(reducedInvariantsVector2);

              // evaluation CC(C_ab,C_cd+h)
              // perturb Cbar
              Cbar = CBarReference;
              Cbar[d][c] += h;
              Cbar[c][d] += h;

              // compute Ibar
              Ibar[0] = Cbar[0][0] + Cbar[1][1] + Cbar[2][2];

              trC2 = 0;
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  trC2 += Cbar[i][j] * Cbar[j][i];
                }
              }
              Ibar[1] = 0.5*((Ibar[0]*Ibar[0]) - trC2);
              Ibar[2] = 1;

              // I4 = a0 • C a0;
              a0Ca0 = 0;
              for (int i = 0; i < 3; i++)
              {
                double ca0_i = 0;
                for (int j = 0; j < 3; j++)
                {
                  ca0_i += Cbar[j][i] * fiberDirection[j];
                }
                a0Ca0 += fiberDirection[i] * ca0_i;
              }

              Ibar[3] = a0Ca0;

              std::vector<double> reducedInvariantsVector3(Ibar.begin(), Ibar.end());
              double Psi01 = PsiExpression.apply(reducedInvariantsVector3);

              double derivative = (Psi11 - Psi01 - Psi10 + PsiReference) / (h*h);
              CCbar2[c][d][b][a] = derivative * 4 * pow(JReference, -4/3.);

              LOG(DEBUG) << a << b << c << d << " Psi: " << Psi11 << " " << Psi01 << " " << Psi10 << " " << PsiReference
                << " -> CCbar2 from Psi: " << derivative << "(ref: " << CCbarReference[c][d][b][a] << ")";

            }  // d
          }  // c
        }  // b
      }  // a

    }

    // compare CCbar2 and CCbarReference
    LOG(DEBUG) << "CCbar:          " << CCbar;
    //LOG(DEBUG) << "CCbar2:         " << CCbar2;
    LOG(DEBUG) << "CCbarReference: " << CCbarReference;

    double error = 0;
    double error2 = 0;

    std::stringstream s,s2,s3;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            s << CCbar[d][c][b][a] << " ";
            s3 << CCbar2[d][c][b][a] << " ";
            s2 << CCbarReference[d][c][b][a] << " ";
            error += MathUtility::sqr(CCbar[d][c][b][a] - CCbarReference[d][c][b][a]);
            error2 += MathUtility::sqr(CCbar2[d][c][b][a] - CCbarReference[d][c][b][a]);
          }  // d
        }  // c
      }  // b
      s << std::endl;
      s2 << std::endl;
      s3 << std::endl;
    }  // a

    LOG(DEBUG) << "CCBar:          " << s.str();
    //LOG(DEBUG) << "CCBar2:         " << s3.str();
    LOG(DEBUG) << "CCbarReference: " << s2.str();

    LOG(DEBUG) << "CCbar error:  " << error/27.;
    //LOG(DEBUG) << "CCbar2 error: " << error2/27.;
    LOG(DEBUG) << "--";

    //if (error > 1)
    //  LOG(FATAL) << "CCbar error is too high";
  }
#endif
#if 0
  if (true)
  {
    double h = 1e-8;

    // compute CCisoReference
    Tensor4<D> CCReference;
    Tensor4<D> CCbarReference;
    Tensor4<3> CCIsoReference;
    computeElasticityTensor(CReference, inverseRightCauchyGreen, JReference, pressure, IbarReference, SBarReference, SIsoReference, fiberDirection,
                            CCbarReference, CCIsoReference, CCReference);

    Tensor2<3> PSbarReference;   //P:Sbar
    PSbarReference = computePSbar(SBarReference, CReference);

    Tensor2<3> SisoReference2 = PSbarReference;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        SisoReference2[b][a] = pow(JReference, -2./3) * PSbarReference[b][a];
      }
    }
/*
    LOG(DEBUG) << "SIsoReference: " << SIsoReference;
    LOG(DEBUG) << "SisoReference2: " << SisoReference2;
    LOG(DEBUG) << "from PSbarReference: " << PSbarReference << ", factor J^-2/3: " << pow(JReference, -2./3);
*/
    double SbarCReference = computeSbarC(SBarReference, CReference);

    // compute CCiso = 2 * ∂Siso/∂C and compare with CCiso
    Tensor4<3> CCiso;
    Tensor4<3> CCiso2;
    Tensor4<3> CCiso3;
    Tensor4<3> T1;
    Tensor4<3> T2;
    Tensor4<3> T2alt1;
    Tensor4<3> A1;
    Tensor4<3> A2;

    // loop over entries of Cbar
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb C
        Tensor2<3> C = CReference;
        C[b][a] += h/2.;
        C[a][b] += h/2.;

        // compute J
        double J2 = MathUtility::computeDeterminant<3>(C);
        double J = sqrt(J2);

        // compute I
        std::array<double,5> I;
        I[0] = C[0][0] + C[1][1] + C[2][2];

        double trC2 = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            trC2 += C[i][j] * C[j][i];
          }
        }
        I[1] = 0.5*((I[0]*I[0]) - trC2);
        I[2] = 1;

        // I4 = a0 • C a0;
        double a0Ca0 = 0;
        for (int i = 0; i < 3; i++)
        {
          double ca0_i = 0;
          for (int j = 0; j < 3; j++)
          {
            ca0_i += C[j][i] * fiberDirection[j];
          }
          a0Ca0 += fiberDirection[i] * ca0_i;
        }

        I[3] = a0Ca0;
        std::array<double,5> Ibar = computeReducedInvariants(I, J);

        // compute C^-1
        double determinant;
        Tensor2<3> CInv = MathUtility::computeInverse<3>(C, determinant);

        // compute Sbar
        Tensor2<3> SBar;
        Tensor2<3> Siso;
        this->computePK2Stress(pressure, C, CInv, Ibar, J, fiberDirection,
                               SBar, Siso);


        //LOG(DEBUG) << "  C: " << C << ", Ibar: " << Ibar[0] << "," << Ibar[1] << ", J: " << J << ", Siso: " << Siso;

        // compute ∂C^-1/∂C:
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double dCdC_abcd = (CInv[d][c] - CInvReference[d][c]) / h;
            double dCdC_abcd_alt1 = -1./2 * (CInvReference[a][c] * CInvReference[b][d] + CInvReference[b][c] * CInvReference[a][d]);
            double dCdC_abcd_alt2 = -CInvReference[a][c] * CInvReference[b][d];
            //LOG(DEBUG) << "++ compare ∂C^-1/∂C: " << dCdC_abcd << "," << dCdC_abcd_alt1 << "," << dCdC_abcd_alt2;
          }
        }


        // compute CCiso
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double derivative = (Siso[d][c] - SIsoReference[d][c]) / h * 2;
            CCiso[b][a][d][c] = derivative;

            LOG(DEBUG) << a << b << c << d << " Siso[d][c]: " << Siso[d][c] << " (diff:" << Siso[d][c]-SIsoReference[d][c]
              << ") -> CCiso: " << derivative << " (ref: " << CCIsoReference[b][a][d][c] << ")";
          }  // d
        }  // c

        // compute T1
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double dJdC_ab = (pow(J,-2./3) - pow(JReference,-2./3)) / h;
            T1[b][a][d][c] = 2*PSbarReference[d][c] * dJdC_ab;

            double T1alternative = -2./3 * SIsoReference[d][c] * CInvReference[b][a];
            LOG(DEBUG) << a << b << c << d << " T1: " << T1[b][a][d][c]<< "=" << 2*PSbarReference[d][c] << "*" << dJdC_ab << ", T1alt=" << T1alternative << "=" << SIsoReference[d][c] << "*" << -2./3 * CInvReference[b][a];
          }
        }

        // compute T2
        Tensor2<3> pSbar = computePSbar(SBar, C);

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double dPSbar_dC = (pSbar[d][c] - PSbarReference[d][c]) / h;
            T2[b][a][d][c] = 2*pow(JReference,-2/3.) * dPSbar_dC;   //T1_cdab


            //LOG(DEBUG) << a << b << c << d << " (" << pSbar[d][c]  << "-" << PSbarReference[d][c] << ")/h=" << dPSbar_dC
            //  << ", JReference: " << JReference << ", T2=" << T2[b][a][d][c] << "=" << 2*pow(JReference,-2/3.) << "*" <<dPSbar_dC;
          }
        }


        // compute CCiso2
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            CCiso2[b][a][d][c] = T1[b][a][d][c] + T2[b][a][d][c];

            //LOG(DEBUG) << a << b << c << d << " T1+T2=CCiso2: " << T1[b][a][d][c] << " + " << T2[b][a][d][c]
            //  << " = " << CCiso2[b][a][d][c] << "(CCiso: " << CCiso[b][a][d][c] << ", ref:" << CCIsoReference[b][a][d][c] << ")";
          }  // d
        }  // c

        // compute CCiso3
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double Siso_cd = pow(J,-2./3)*pSbar[d][c];
            //LOG(DEBUG) << "Siso: " << Siso[d][c] << ", " << Siso_cd;

            double value = (Siso_cd - pow(JReference,-2./3)*PSbarReference[d][c]) / h * 2;
            CCiso3[b][a][d][c] = value;

            //LOG(DEBUG) << a << b << c << d << " CCiso3 value: " << value
            //  << " (CCiso: " << CCiso[b][a][d][c] << ", CCiso2: " << CCiso2[b][a][d][c] << ", ref: " << CCIsoReference[b][a][d][c] << ")";
          }  // d
        }  // c


        // compute A1

        // perturb Cbar
        Tensor2<3> Cbar = CBarReference;
        Cbar[b][a] += h/2.;
        Cbar[a][b] += h/2.;

        // compute Ibar
        Ibar[0] = Cbar[0][0] + Cbar[1][1] + Cbar[2][2];

        trC2 = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            trC2 += Cbar[i][j] * Cbar[j][i];
          }
        }
        Ibar[1] = 0.5*((Ibar[0]*Ibar[0]) - trC2);
        Ibar[2] = 1;

        // I4 = a0 • C a0;
        a0Ca0 = 0;
        for (int i = 0; i < 3; i++)
        {
          double ca0_i = 0;
          for (int j = 0; j < 3; j++)
          {
            ca0_i += Cbar[j][i] * fiberDirection[j];
          }
          a0Ca0 += fiberDirection[i] * ca0_i;
        }

        Ibar[3] = a0Ca0;

        // compute Cbar^-1
        Tensor2<3> CbarInv = MathUtility::computeInverse<3>(Cbar, determinant);

        // compute C
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            C[d][c] = pow(JReference,2./3) * Cbar[d][c];
          }
        }
        CInv = MathUtility::computeInverse<3>(C, determinant);


        // compute Sbar
        this->computePK2Stress(pressure, C, CInv, Ibar, JReference, fiberDirection,
                                                        SBar, Siso);

        //LOG(DEBUG) << "  Cbar: " << Cbar << ", Ibar: " << Ibar[0] << "," << Ibar[1] << ", SBar: " << SBar;

        // compute ∂Cbar^-1/∂Cbar and ∂C^-1/∂Cbar

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double dCbardC_abcd = (CbarInv[d][c] - CBarInvReference[d][c]) / h;
            double dCbardC_abcd_alt1 = -1./2 * (CBarInvReference[a][c] * CBarInvReference[b][d] + CBarInvReference[b][c] * CBarInvReference[a][d]);
            double dCbardC_abcd_alt2 = -CBarInvReference[a][c] * CBarInvReference[b][d];
            //LOG(DEBUG) << "compare ∂Cbar^-1/∂Cbar: " << dCbardC_abcd << "," << dCbardC_abcd_alt1 << "," << dCbardC_abcd_alt2;

            double dCdC_abcd = (CInv[d][c] - CInvReference[d][c]) / h;
            double dCdC_abcd_alt1 = dCbardC_abcd_alt1 * pow(JReference,-2./3);
            double dCdC_abcd_alt2 = dCbardC_abcd_alt2 * pow(JReference,-2./3);
            //LOG(DEBUG) << "compare ∂C^-1/∂Cbar: " << dCdC_abcd << "," << dCdC_abcd_alt1 << "," << dCdC_abcd_alt2;
          }
        }

        // compute CCbar = 2J^{-4/3} * ∂Sbar / ∂Cbar
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double CCbar1_cdab = (SBar[d][c] - SBarReference[d][c]) / h * 2 * pow(JReference, -4./3);
            double CCbarRef = CCbarReference[b][a][d][c];

            std::string s("");
            if (fabs(CCbar1_cdab - CCbarRef) > 1e-4)
              s = " *";

            LOG(DEBUG) << c << d << a << b << ": CCbar1: "<< CCbar1_cdab << ", CCbarRef: " << CCbarRef << s;
          }
        }

        // compute A1
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double derivative = (computeSbarC(SBar, C) - SbarCReference) / h;

            A1[b][a][d][c] = 2*pow(JReference,-4/3.) * CInvReference[d][c] * derivative;   //A1_cdab

            // A1alt
            double cCcbar_ab = 0;
            for (int e = 0; e < 3; e++)
            {
              for (int f = 0; f < 3; f++)
              {
                cCcbar_ab += CReference[f][e] * CCbarReference[b][a][f][e];   // C_ef * CC_efab
              }
            }
            double A1alt_cdab = CInvReference[d][c] * (cCcbar_ab + 2*pow(JReference,-2./3)*SBarReference[b][a]);
            //LOG(DEBUG) << a << b << c << d << " A1: " << A1[b][a][d][c] << ", alt: " << A1alt_cdab;
          }
        }

        // compute A2
        double det;
        CInv = MathUtility::computeInverse<3>(C, det);
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double derivative = (CInv[d][c] - CInvReference[d][c]) / h;

            A2[b][a][d][c] = 2*pow(JReference,-4/3.) * SbarCReference * derivative;   //A2_cdab

            // A2alt
            // cdab -> ca db  cb da
            double c1c1_cdab = 1./2*
              (CInvReference[a][c]*CInvReference[b][d] + CInvReference[b][c]*CInvReference[a][d]);

            //c1c1_cdab = (CInvReference[a][c]*CInvReference[b][d]);

              double A2alt_cdab = -2*pow(JReference,-2./3) * SbarCReference * c1c1_cdab;
            LOG(DEBUG) << a << b << c << d << " A2: " << A2[b][a][d][c] << ", alt: " << A2alt_cdab;
          }
        }

      }  // b
    }  // a


    // compute T2alt1

    // loop over entries of T2al1
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int e = 0; e < 3; e++)
        {
          for (int f = 0; f < 3; f++)
          {
            T2alt1[f][e][b][a] = 0.0;
            for (int c = 0; c < 3; c++)
            {
              for (int d = 0; d < 3; d++)
              {
                double leftTerm_abcd = (CCbarReference[d][c][b][a] - 1./3 * (A1[d][c][b][a] + A2[d][c][b][a]));

                int delta_ce = (c == e? 1 : 0);
                int delta_df = (d == f? 1 : 0);
                const int Ii = delta_ce * delta_df;       // II = (δ_AC*δ_BD + δ_AD*δ_BC) / 2

                const double Cc = CInvReference[f][e] * CReference[d][c];     // CC = C^{-1}_AB * C_CD
                const double pp_efcd = (Ii - 1./3 * Cc);

                double rightTerm_cdef = pp_efcd;

                T2alt1[f][e][b][a] += leftTerm_abcd * rightTerm_cdef;
              }
            }
          }
        }

      }  // b
    }  // a

    // compute T2alt2
    Tensor4<3> T2alt2;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            // Term 1
            double term1_abcd = 0;
            for (int g = 0; g < 3; g++)
            {
              for (int h = 0; h < 3; h++)
              {
                double iccCbar_abgh = 0;
                for (int e = 0; e < 3; e++)
                {
                  for (int f = 0; f < 3; f++)
                  {
                    int delta_ae = (a == e? 1 : 0);
                    int delta_bf = (b == f? 1 : 0);
                    double icc_abef = delta_ae * delta_bf - 1./3 * CInvReference[b][a] * CReference[f][e];

                    iccCbar_abgh += icc_abef * CCbarReference[h][g][f][e];
                  }
                }

                int delta_cg = (c == g? 1 : 0);
                int delta_dh = (d == h? 1 : 0);
                const int Ii = delta_cg * delta_dh;       // II = (δ_AC*δ_BD + δ_AD*δ_BC) / 2

                const double Cc = CInvReference[d][c] * CReference[h][g];     // CC = C^{-1}_AB * C_CD
                const double pp_cdgh = (Ii - 1./3 * Cc);

                double rightTerm_ghcd = pp_cdgh;    // transpose

                term1_abcd += iccCbar_abgh * rightTerm_ghcd;
              }  // h
            }  // g

            // Term 2
            double SBarP_cd = 0;
            for (int e = 0; e < 3; e++)
            {
              for (int f = 0; f < 3; f++)
              {
                int delta_ce = (c == e? 1 : 0);
                int delta_df = (d == f? 1 : 0);
                const int Ii = delta_ce * delta_df;       // II = (δ_AC*δ_BD + δ_AD*δ_BC) / 2

                const double Cc = CInvReference[d][c] * CReference[f][e];     // CC = C^{-1}_AB * C_CD
                const double pp_cdef = (Ii - 1./3 * Cc);

                double rightTerm_efcd = pp_cdef;    // transpose

                SBarP_cd += SBarReference[f][e] * rightTerm_efcd;
              }
            }

            //LOG(DEBUG) << "term2 : " << pow(JReference,-2./3) * SBarP_cd << "," << SIsoReference[d][c];

            double term2_abcd = -2/3.*CInvReference[b][a] * pow(JReference,-2./3) * SBarP_cd;

            // Term 3
            double term3_prefactor = 0;
            for (int e = 0; e < 3; e++)
            {
              for (int f = 0; f < 3; f++)
              {
                term3_prefactor += 2./3 * pow(JReference, -2./3) * SBarReference[f][e] * CReference[f][e];
              }
            }

            double term3_abcd = 0;
            for (int e = 0; e < 3; e++)
            {
              for (int f = 0; f < 3; f++)
              {
                double cc_abef = 1./2 *
                  (CInvReference[e][a] * CInvReference[f][b] + CInvReference[f][a]*CInvReference[e][b]);
                //cc_abef = CInvReference[e][a] * CInvReference[f][b];

                int delta_ce = (c == e? 1 : 0);
                int delta_df = (d == f? 1 : 0);
                const int Ii = delta_ce * delta_df;       // II = (δ_AC*δ_BD + δ_AD*δ_BC) / 2

                const double Cc = CInvReference[d][c] * CReference[f][e];     // CC = C^{-1}_AB * C_CD
                const double pp_cdef = (Ii - 1./3 * Cc);

                double rightTerm_efcd = pp_cdef;    // transpose

                term3_abcd += cc_abef * rightTerm_efcd;
              }
            }


            double term3alt_abcd = CInvReference[c][a] * CInvReference[d][b] - 1./3*CInvReference[b][a]*CInvReference[d][c];
            //LOG(DEBUG) << "term 3 : " << term3_abcd << "," << term3alt_abcd;

            term3_abcd *= term3_prefactor;


            T2alt2[d][c][b][a] = term1_abcd + term2_abcd + term3_abcd;
          }  // d
        }  // c
      }   // b
    }  // a



    // compare T2
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            LOG(DEBUG) << "compare T2: T2_" << a << b << c << d << ": " << T2[d][c][b][a] << "," << T2alt1[d][c][b][a] << "," << T2alt2[d][c][b][a];
          }
        }
      }   // b
    }  // a

    // compute CCiso4
    Tensor4<3> CCiso4;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            // 1.Term P:Cbar:P^T
            double term1_abcd = 0;
            for (int g = 0; g < 3; g++)
            {
              for (int h = 0; h < 3; h++)
              {
                double pCbar_abgh = 0;
                for (int e = 0; e < 3; e++)
                {
                  for (int f = 0; f < 3; f++)
                  {
                    int delta_ae = (a == e? 1 : 0);
                    int delta_bf = (b == f? 1 : 0);
                    double icc_abef = delta_ae * delta_bf - 1./3 * CInvReference[b][a] * CReference[f][e];

                    pCbar_abgh += icc_abef * CCbarReference[h][g][f][e];
                  }
                }

                int delta_cg = (c == g? 1 : 0);
                int delta_dh = (d == h? 1 : 0);
                const int Ii = delta_cg * delta_dh;       // II = (δ_AC*δ_BD + δ_AD*δ_BC) / 2

                const double Cc = CInvReference[d][c] * CReference[h][g];     // CC = C^{-1}_AB * C_CD
                const double pp_cdgh = (Ii - 1./3 * Cc);

                double rightTerm_ghcd = pp_cdgh;    // transpose

                term1_abcd += pCbar_abgh * rightTerm_ghcd;
              }  // h
            }  // g

            // 2.Term 2/3 Tr(J^-2/3 Sbar) Ptilde
            double TrSbar = 0;
            for (int e = 0; e < 3; e++)
            {
              for (int f = 0; f < 3; f++)
              {
                TrSbar += SBarReference[f][e] * CReference[f][e];
              }
            }

            double cInvCInv = 0.5*(CInvReference[c][a] * CInvReference[d][b] + CInvReference[d][a] * CInvReference[c][b]);
            double pTilde_abcd = cInvCInv - 1./3*CInvReference[b][a]*CInvReference[d][c];
            double term2_abcd = 2./3 * pow(JReference, -2./3) * TrSbar * pTilde_abcd;

            // 3.Term -2/3 (C^-1 dyad Siso + Siso dyad C^-1)
            double term3_abcd = -2/3.*(CInvReference[b][a] * SIsoReference[d][c] + SIsoReference[b][a] * CInvReference[d][c]);

            CCiso4[d][c][b][a] = term1_abcd + term2_abcd + term3_abcd;

            //LOG(DEBUG) << a << b << c << d << " CCiso correct: " << CCiso4[d][c][b][a] << ": " << term1_abcd << "," << term2_abcd << "," << term3_abcd;
            LOG(DEBUG) << a << b << c << d << " CCiso4 value: " << CCiso4[d][c][b][a]
              << " (CCiso: " << CCiso[d][c][b][a] << ", CCiso2: " << CCiso2[d][c][b][a] << ", ref: " << CCIsoReference[d][c][b][a] << ")";
          }  // d
        }  // c
      }   // b
    }  // a


    // compare CCiso and CCIsoReference
    //LOG(DEBUG) << "CCiso:          " << CCiso;
    //LOG(DEBUG) << "CCiso2:         " << CCiso2;
    //LOG(DEBUG) << "CCIsoReference: " << CCIsoReference;

    double error = 0;
    double error2 = 0;
    double error3 = 0;

    std::stringstream s,s2,s3,s4;
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            s << CCiso[d][c][b][a] << " ";
            s3 << CCiso2[d][c][b][a] << " ";
            s4 << CCiso4[d][c][b][a] << " ";
            s2 << CCIsoReference[d][c][b][a] << " ";
            error += MathUtility::sqr(CCiso[d][c][b][a] - CCIsoReference[d][c][b][a]);
            error2 += MathUtility::sqr(CCiso2[d][c][b][a] - CCIsoReference[d][c][b][a]);
            error3 += MathUtility::sqr(CCiso4[d][c][b][a] - CCIsoReference[d][c][b][a]);
          }  // d
        }  // c
      }  // b
      s << std::endl;
      s2 << std::endl;
      s3 << std::endl;
      s4 << std::endl;
    }  // a

    LOG(DEBUG) << "CCiso:          " << s.str();
    LOG(DEBUG) << "CCiso2:         " << s3.str();
    LOG(DEBUG) << "CCiso4:         " << s4.str();
    LOG(DEBUG) << "CCIsoReference: " << s2.str();

    LOG(DEBUG) << "CCiso error:  " << error/27.;
    LOG(DEBUG) << "CCiso2 error: " << error2/27.;
    LOG(DEBUG) << "CCiso4 error: " << error3/27.;
    LOG(DEBUG) << "--";

  }
#endif

  if (true)
  {
    double h = 1e-6;

    // compute CCisoReference
    Tensor4<3> CCReference;
    Tensor4<3> CCbarReference;
    Tensor4<3> CCIsoReference;
    computeElasticityTensor(CReference, CInvReference, JReference, pressure, IbarReference, SBarReference, SIsoReference, fiberDirection,
                            CCbarReference, CCIsoReference, CCReference);

    bool match = true;

    // loop over entries of C
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb C
        Tensor2<3> C = CReference;
        C[b][a] += h/2.;
        C[a][b] += h/2.;

        // compute J
        double J2 = MathUtility::computeDeterminant<3>(C);
        double J = sqrt(J2);

        // compute I
        std::array<double,5> I;
        I[0] = C[0][0] + C[1][1] + C[2][2];

        double trC2 = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            trC2 += C[i][j] * C[j][i];
          }
        }
        I[1] = 0.5*((I[0]*I[0]) - trC2);
        I[2] = 1;

        // I4 = a0 • C a0;
        double a0Ca0 = 0;
        for (int i = 0; i < 3; i++)
        {
          double ca0_i = 0;
          for (int j = 0; j < 3; j++)
          {
            ca0_i += C[j][i] * fiberDirection[j];
          }
          a0Ca0 += fiberDirection[i] * ca0_i;
        }

        I[3] = a0Ca0;
        std::array<double,5> Ibar = computeReducedInvariants(I, J);

        // compute C^-1
        double determinant;
        Tensor2<3> CInv = MathUtility::computeInverse<3>(C, determinant);

        // compute Sbar and Siso
        Tensor2<3> SBar;
        Tensor2<3> Siso;
        this->computePK2Stress(pressure, C, CInv, Ibar, J, fiberDirection,
                               SBar, Siso);

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {
            double derivative = (Siso[d][c] - SIsoReference[d][c]) / h;
            double CCiso_cdab = derivative * 2;

            std::string s("");
            if (fabs(CCiso_cdab - CCIsoReference[b][a][d][c]) > 1e-3)
            {
              match = false;
              double factor = CCiso_cdab / CCIsoReference[b][a][d][c];
              std::stringstream ss;
              ss << "* (" << factor << ")";
              s = ss.str();
            }

            //LOG(DEBUG) << c << d << a << b << ": CCiso numeric: " << CCiso_cdab << ", analytic: " << CCIsoReference[b][a][d][c] << " " << s;
          }  // d
        }  // c

      }  // b
    }  // a

    LOG(DEBUG) << "match: " << match;
  }

#endif
}

} // namespace
