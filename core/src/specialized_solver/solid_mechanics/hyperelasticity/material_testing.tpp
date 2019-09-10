#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "equation/mooney_rivlin_incompressible.h"

namespace SpatialDiscretization
{

template<typename Term>
void HyperelasticitySolver<Term>::
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
  LOG(DEBUG) << "materialTesting, parameters: " << this->materialParameters_ << ", C: " << rightCauchyGreen;

  // set reference values
  double JReference                  = deformationGradientDeterminant;
  Tensor2<3> CReference              = rightCauchyGreen;
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

  // compute Sbar from input values
  Tensor2<3> SBarReference;
  Tensor2<3> Siso;
  this->computePK2Stress(pressure, CReference, inverseRightCauchyGreen, IbarReference, JReference, fiberDirection,
                                                  SBarReference, Siso
                                                );

  // compute Psi
  auto PsiExpression = Term::strainEnergyDensityFunctionIsochoric;
  std::vector<double> reducedInvariantsVector(IbarReference.begin(), IbarReference.end());
  const double PsiReference = PsiExpression.apply(reducedInvariantsVector);


  LOG(DEBUG) << "reference values: ";
  LOG(DEBUG) << "  J: " << JReference << ", C: " << CReference;
  LOG(DEBUG) << "  Ibar: " << IbarReference[0] << "," << IbarReference[1] << ", CBar: " << CBarReference;
  LOG(DEBUG) << "  Psi: " << PsiReference;
  LOG(DEBUG) << "  Sbar: " << SBarReference;

  // test if Sbar = 2 ∂ψ/∂Cbar by evaluating the derivative numerically (finite differences)
  Tensor2<3> SPsi;

  double h = 1e-10;

  if (false)
  {
    // loop over entries of Cbar
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb Cbar
        Tensor2<3> Cbar = CBarReference;
        Cbar[b][a] += h;

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
        LOG(DEBUG) << "  Cbar: " << Cbar << ", Ibar: " << Ibar[0] << "," << Ibar[1] << ", Psi: " << Psi;
      }  // b
    }  // a
  }

  // compare SPsi and SBarReference
  if (false)
  {
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

    if (error > 1)
      LOG(ERROR) << "Sbar error is too high";
  }

  // -----------------------------------

  const int D = 3;
  if (false)
  {
    // compute reference C
    Tensor4<D> CC;
    Tensor4<D> CCbarReference;
    computeElasticityTensor(CReference, inverseRightCauchyGreen, JReference, pressure, IbarReference, SBarReference, Siso, fiberDirection,
                            CCbarReference, CC);

    Tensor4<D> CCbar;

    // test if Cbar = J^{-4/3} * 2 * ∂Sbar/∂Cbar by evaluating the 1st derivative numerically
    // loop over entries of Cbar
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb Cbar
        Tensor2<3> Cbar = CBarReference;
        Cbar[b][a] += h;

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
        double determinant;
        Tensor2<3> CbarInv = MathUtility::computeInverse<3>(Cbar, determinant);

        // compute Sbar, using J=1 such that input C = Cbar
        Tensor2<3> SBar;
        Tensor2<3> Siso;
        this->computePK2Stress(pressure, Cbar, CbarInv, Ibar, 1, fiberDirection,
                                                        SBar, Siso);

        LOG(DEBUG) << "  Cbar: " << Cbar << ", Ibar: " << Ibar[0] << "," << Ibar[1] << ", SBar: " << SBar;

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {

            double derivative = (SBar[d][c] - SBarReference[d][c]) / h * 2 * pow(JReference, -4/3.);
            CCbar[b][a][d][c] = derivative;

            LOG(DEBUG) << a << b << c << d << " SBar[d][c]: " << SBar[d][c] << " (diff:" << SBar[d][c]-SBarReference[d][c]
              << ") -> value: " << derivative << "(ref: " << CCbarReference[d][c][b][a] << ")";
          }  // d
        }  // c
      }  // b
    }  // a


    // test if Cbar = J^{-4/3} * 4 * ∂^2Sbar/∂Cbar^2 by evaluating the 2nd derivative numerically

    Tensor4<D> CCbar2;
    if (true)
    {
      h = 1e-3;

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
              Cbar[d][c] += h;

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
                << " -> -> value: " << derivative << "(ref: " << CCbarReference[c][d][b][a] << ")";

            }  // d
          }  // c
        }  // b
      }  // a

    }

    // compare CCbar2 and CCbarReference
    LOG(DEBUG) << "CCbar:          " << CCbar;
    LOG(DEBUG) << "CCbar2:         " << CCbar2;
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
    LOG(DEBUG) << "CCBar2:          " << s3.str();
    LOG(DEBUG) << "CCbarReference: " << s2.str();

    LOG(DEBUG) << "CCbar error: " << error/27.;
    LOG(DEBUG) << "CCbar2 error: " << error2/27.;
    LOG(DEBUG) << "--";

    if (error > 1)
      LOG(FATAL) << "CCbar error is too high";
  }

  if (true)
  {
    // compute CCisoReference
    Tensor4<D> CC;
    Tensor4<D> CCbarReference;
    computeElasticityTensor(CReference, inverseRightCauchyGreen, JReference, pressure, IbarReference, SBarReference, Siso, fiberDirection,
                            CCbarReference, CC);Todo: return CCisoReference

    // compute CCiso = 2 * ∂Siso/∂C and compare with CCiso
    Tensor2<3> CCiso;
    // loop over entries of Cbar
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        // perturb Cbar
        Tensor2<3> C = CReference;
        C[b][a] += h;

        // compute J
        double J2 = MathUtility::computeDeterminant<3>(C);
        double J = sqrt(J2);

        // compute I
        std::array<double,5> I
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

        // compute C^-1
        double determinant;
        Tensor2<3> CInv = MathUtility::computeInverse<3>(C, determinant);

        // compute Sbar, using J=1 such that input C = Cbar
        Tensor2<3> SBar;
        Tensor2<3> Siso;
        this->computePK2Stress(pressure, C, CInv, I, J, fiberDirection,
                                                        SBar, Siso);

        LOG(DEBUG) << "  C: " << C << ", I: " << I[0] << "," << I[1] << ", Siso: " << Siso;

        for (int c = 0; c < 3; c++)
        {
          for (int d = 0; d < 3; d++)
          {

            double derivative = (Siso[d][c] - SisoReference[d][c]) / h * 2;
            CCiso[b][a][d][c] = derivative;

            LOG(DEBUG) << a << b << c << d << " Siso[d][c]: " << Siso[d][c] << " (diff:" << Siso[d][c]-SisoReference[d][c]
              << ") -> value: " << derivative << "(ref: " << CCisoReference[b][a][d][c] << ")";
          }  // d
        }  // c
      }  // b
    }  // a
  }
}

} // namespace
