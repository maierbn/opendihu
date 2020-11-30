#pragma once

#include "equation/static.h"

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace Equation
{
namespace SolidMechanics
{

/** Base class that defines the parameters of the strain energy density function.
 *  All hyperelasticity formulations should inherit from this class.
 */
struct HyperelasticityBase
{
  // define helper variables for strain energy density function
  // note, the order of VAR(0) - VAR(16) is reflected in the code and cannot be changed anymore

  // reduced invariants, arguments of `strainEnergyDensityFunctionIsochoric`
  static constexpr auto Ibar1 = VAR(0); //< 1st reduced or modified strain invariant Ibar1 = tr(Cbar) = J^{-2/3}*I_1
  static constexpr auto Ibar2 = VAR(1); //< 2nd reduced or modified strain invariant Ibar2 = 1/2 (tr(Cbar)^2 - tr(Cbar^2)) = J^{-2/3}*I_2
                                        // Note, there is no 3rd reduced or modified strain invariant needed, Ibar3 = det(Cbar) = 1 (incompressibility)
  static constexpr auto Ibar4 = VAR(2); //< 4th reduced or modified strain invariant Ibar4 = a0•Cbar a0
  static constexpr auto Ibar5 = VAR(3); //< 5th reduced or modified strain invariant Ibar5 = a0•Cbar^2 a0
  static constexpr auto lambda = sqrt(Ibar4);    //< fiber stretch, helper variable that can also be used in `strainEnergyDensityFunctionIsochoric`

  // volume factor, argument of `strainEnergyDensityFunctionVolumetric` (only for compressible material)
  static constexpr auto J = VAR(4);     //< volume factor, J = det(F), only for compressible material (otherwise it is 1)

  // invariants, arguments of `strainEnergyDensityFunctionCoupled`
  static constexpr auto I1 = VAR(5);    //< 1st strain invariant I1 = tr(C)
  static constexpr auto I2 = VAR(6);    //< 2nd strain invariant I2 = 1/2 (tr(C)^2 - tr(C^2))
  static constexpr auto I3 = VAR(7);    //< 3rd strain invariant I3 = det(C) = J^2

  // components of the right Cauchy Green tensor, arguments of `strainEnergyDensityFunctionCoupledDependentOnC
  static constexpr auto C11 = VAR(8);   //< entry C11 of the right Cauchy Green tensor, C
  static constexpr auto C12 = VAR(9);   //< entry C12 = C21 of the right Cauchy Green tensor, C
  static constexpr auto C13 = VAR(10);  //< entry C13 = C31 of the right Cauchy Green tensor, C
  static constexpr auto C22 = VAR(11);  //< entry C22 of the right Cauchy Green tensor, C
  static constexpr auto C23 = VAR(12);  //< entry C23 = C32 of the right Cauchy Green tensor, C
  static constexpr auto C33 = VAR(13);  //< entry C33 of the right Cauchy Green tensor, C

  static constexpr auto a1 = VAR(14);   //< entry a0_1 of the fiber direction, a0
  static constexpr auto a2 = VAR(15);   //< entry a0_2 of the fiber direction, a0
  static constexpr auto a3 = VAR(16);   //< entry a0_3 of the fiber direction, a0

  static constexpr auto I4
    = C11*a1*a1 + INT(2)*C12*a1*a2
      + INT(2)*C13*a1*a3 + C22*a2*a2
      + INT(2)*C23*a2*a3 + C33*a3*a3;   //< non-reduced 4th strain invariant, I4 = a0•C a0
};

/** Isotropic hyperelastic solid mechanics formulation with arbitrary strain energy density function
 */
struct MooneyRivlinIncompressible3D : HyperelasticityBase
{
  static constexpr bool isIncompressible = true;       //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = false;    //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress

  // material parameters
  static constexpr auto c1 = PARAM(0);   //< material parameter
  static constexpr auto c2 = PARAM(1);   //< material parameter

  static constexpr int nMaterialParameters = 2;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3));
    //= c1*(pow(Ibar1,INT(3)) - INT(3)) + c2*(pow(Ibar2,INT(2))*Ibar1 - INT(3));  // complicated artifical equation for debugging

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

  //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C
  //! it must only depend on variables C11, C12, C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);
};

/** Transversely-isotropic hyperelastic solid mechanics formulation, here Mooney-Rivlin
 */
struct TransverselyIsotropicMooneyRivlinIncompressible3D : HyperelasticityBase
{
  static constexpr bool isIncompressible = true;       //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = true;     //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress

  // material parameters
  static constexpr auto c1 = PARAM(0);     //< material parameter
  static constexpr auto c2 = PARAM(1);     //< material parameter
  static constexpr auto b = PARAM(2);      //< material parameter
  static constexpr auto d = PARAM(3);      //< material parameter

  static constexpr int nMaterialParameters = 4;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  //! function from chemo-electro-mechanical muscle model (Heidlauf 2013 "Modeling the chemo-electro-mechanical behaviour of... p.4)
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3)) + b/d * (pow(lambda, d) - INT(1)) - b*ln(lambda);

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

  //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C
  //! it must only depend on variables C11, C12, C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);
};

//! this is not yet fully implemented
/** Model used in chemo-electro-mechanical model
 */
struct TransverselyIsotropicMooneyRivlinIncompressibleActive3D : HyperelasticityBase
{
  static constexpr bool isIncompressible = true;       //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = true;     //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress = true;       //< if the value of an active stress term will be added to the stress

  // material parameters
  static constexpr auto c1 = PARAM(0);     //< material parameter
  static constexpr auto c2 = PARAM(1);     //< material parameter
  static constexpr auto b = PARAM(2);      //< material parameter
  static constexpr auto d = PARAM(3);      //< material parameter

  static constexpr int nMaterialParameters = 4;  //< number of material parameters

  // define helper variables for strain energy density function
  static constexpr auto lambda = sqrt(Ibar4);    //< fiber stretch

  //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  //! function from chemo-electro-mechanical muscle model (Heidlauf 2013 "Modeling the chemo-electro-mechanical behaviour of... p.4)
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3)) + b/d * (pow(lambda, d) - INT(1)) - b*ln(lambda);

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

  //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C
  //! it must only depend on variables C11, C12, C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);
};

/** Tendon model from "Carniel, T. A., & Fancello, E. A. (2017). A transversely isotropic coupled hyperelastic model for the mechanical behavior of tendons. Journal of biomechanics, 54, 49-57."
 */
struct HyperelasticTendon : HyperelasticityBase
{
  static constexpr bool isIncompressible = false;      //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = true;     //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress

  // material parameters that will be set by the python settings
  static constexpr auto c   = PARAM(0);     //< material parameter
  static constexpr auto ca  = PARAM(1);     //< material parameter
  static constexpr auto ct  = PARAM(2);      //< material parameter
  static constexpr auto cat = PARAM(3);      //< material parameter
  static constexpr auto ctt = PARAM(4);      //< material parameter
  static constexpr auto mu  = PARAM(5);      //< material parameter
  static constexpr auto k1  = PARAM(6);      //< material parameter
  static constexpr auto k2  = PARAM(7);      //< material parameter

  static constexpr int nMaterialParameters = 8;  //< number of material parameters

  static constexpr auto psi_NH = mu/INT(2) * (Ibar1 - INT(3));     //< Neo-Hookean expression
  //static constexpr auto psi_NH = INT(0);

  static constexpr auto E11 = INT(1)/INT(2) * ln(C11);        //< components of the logarithmic strain measure, E_(0) = ln(U) = 1/2*ln(C)
  static constexpr auto E12 = INT(1)/INT(2) * ln(C12);
  static constexpr auto E13 = INT(1)/INT(2) * ln(C13);
  static constexpr auto E22 = INT(1)/INT(2) * ln(C22);
  static constexpr auto E23 = INT(1)/INT(2) * ln(C23);
  static constexpr auto E33 = INT(1)/INT(2) * ln(C33);
  static constexpr auto Q = ca*E11*E11
                            + ct*E22*E22 + ct*E33*E33
                            + cat*E11*E22*INT(2) + cat*E11*E33*INT(2)
                            + ctt*E22*E33*INT(2);                  //< exponent, Q = E:A:E with A (2nd order tensor) given by ca,ct,cat,ctt
  static constexpr auto psi_Fung = c/INT(2) * (exp(Q) - INT(1));   //< coupled orthotropic model of Fung at al. (1979)
//  static constexpr auto psi_Fung = INT(0);

  static constexpr auto psi_f = k1 * pow(I4 - INT(1), INT(2)) + k2 * pow(I4 - INT(1), INT(3));    //< fiber strain energy

  //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric = psi_NH;

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

  //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C
  //! it must only depend on variables C11, C12, C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = psi_Fung + psi_f;
};

}  // namespace
}  // namespace
