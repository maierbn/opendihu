#pragma once

#include "equation/static.h"

#include "semt/Semt.h"
#include "semt/Shortcuts.h"
#include "equation/mooney_rivlin_incompressible.h"

namespace Equation
{
namespace SolidMechanics
{

/** Saint Venant-Kirchhoff model Psi = 1/2 λ (trE)^2 + μ E^2
 */
struct SaintVenantKirchhoff : HyperelasticityBase
{
  static constexpr bool isIncompressible = false;       //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = false;    //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress

  // material parameters
  static constexpr auto lambd = PARAM(0);   //< material parameter
  static constexpr auto mu = PARAM(1);      //< material parameter

  static constexpr int nMaterialParameters = 2;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric = INT(0);

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  // I1 = tr(C), tr(E) = 1/2(tr(C) - tr(I)) = 1/2*(tr(C) - 3) = 1/2*(I1-3)
  static constexpr auto trE = INT(1)/INT(2) * (I1 - INT(3));

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(1)/INT(2) * lambd * pow(trE,INT(2));

  // compute components of E = 1/2(C - I)
  static constexpr auto E11 = INT(1)/INT(2) * (C11 - INT(1));
  static constexpr auto E12 = INT(1)/INT(2) * C12;
  static constexpr auto E13 = INT(1)/INT(2) * C13;
  static constexpr auto E22 = INT(1)/INT(2) * (C22 - INT(1));
  static constexpr auto E23 = INT(1)/INT(2) * C23;
  static constexpr auto E33 = INT(1)/INT(2) * (C33 - INT(1));

  // compute E^2 = E:E
  static constexpr auto EE = pow(E11,INT(2)) + INT(2)*pow(E12,INT(2)) + INT(2)*pow(E13,INT(2)) + pow(E22,INT(2)) + INT(2)*pow(E23,INT(2)) + pow(E33,INT(2));

  //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C
  //! it must only depend on variables C11, C12, C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = mu * EE;
};


}  // namespace
}  // namespace
