#pragma once

#include "equation/static.h"

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace Equation
{
namespace SolidMechanics
{

/** Isotropic hyperelastic solid mechanics formulation with arbitrary strain energy density function
 */
struct MooneyRivlinIncompressible3D
{
  static constexpr bool usesFiberDirection = false;   //< if the equation depends on a fiber direction, i.e. has 4th and 5th invariant
  static constexpr bool usesActiveStress = false;   //< if the equation adds an active stress value to the stress

  // define helper variables for strain energy density function
  // reduced invariants
  static constexpr auto Ibar1 = VAR(0);    //< 1st reduced invariant
  static constexpr auto Ibar2 = VAR(1);    //< 2nd reduced invariant
  static constexpr auto Ibar4 = VAR(4);    //< 4th reduced invariant (not used)
  static constexpr auto J = VAR(2);        //< jacobian J=det F

  // parameters
  static constexpr auto c1 = PARAM(0);   //< material parameter
  static constexpr auto c2 = PARAM(1);   //< material parameter

  static constexpr int nMaterialParameters = 2;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy density function, Psi_iso, in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3));
    //= c1*(pow(Ibar1,INT(3)) - INT(3)) + c2*(pow(Ibar2,INT(2))*Ibar1 - INT(3));  // complicated artifical equation for debugging
};

/** Transversely-isotropic hyperelastic solid mechanics formulation, here Mooney-Rivlin
 */
struct TransverselyIsotropicMooneyRivlinIncompressible3D
{
  static constexpr bool usesFiberDirection = true;   //< if the equation depends on a fiber direction, i.e. has 4th and 5th invariant
  static constexpr bool usesActiveStress = false;   //< if the equation adds an active stress value to the stress

  // define helper variables for strain energy density function
  // reduced invariants
  static constexpr auto Ibar1 = VAR(0);    //< 1st reduced invariant
  static constexpr auto Ibar2 = VAR(1);    //< 2nd reduced invariant
  static constexpr auto Ibar3 = VAR(2);    //< 3rd reduced invariant
  static constexpr auto Ibar4 = VAR(3);    //< 4th reduced invariant
  static constexpr auto lambda = sqrt(Ibar4);    //< fiber stretch
  static constexpr auto J = VAR(2);        //< jacobian J=det F

  // parameters
  static constexpr auto c1 = PARAM(0);   //< material parameter
  static constexpr auto c2 = PARAM(1);   //< material parameter
  static constexpr auto b = PARAM(2);   //< material parameter
  static constexpr auto d = PARAM(3);   //< material parameter

  static constexpr int nMaterialParameters = 4;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy density function, Psi_iso, in terms of the reduced invariants
  // function from chemo-electro-mechanical muscle model (Heidlauf 2013 "Modeling the chemo-electro-mechanical behaviour of... p.4)
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3)) + b/d * (pow(lambda, d) - INT(1)) - b*ln(lambda);
};

//! this is not yet fully implemented
/** Model used in chemo-electro-mechanical model
 */
struct TransverselyIsotropicMooneyRivlinIncompressibleActive3D
{
  static constexpr bool usesFiberDirection = true;   //< if the equation depends on a fiber direction, i.e. has 4th and 5th invariant
  static constexpr bool usesActiveStress = true;   //< if the equation adds an active stress value to the stress

  // define helper variables for strain energy density function
  // reduced invariants
  static constexpr auto Ibar1 = VAR(0);    //< 1st reduced invariant
  static constexpr auto Ibar2 = VAR(1);    //< 2nd reduced invariant
  static constexpr auto Ibar4 = VAR(4);    //< 4th reduced invariant
  static constexpr auto lambda = sqrt(Ibar4);    //< fiber stretch
  static constexpr auto J = VAR(2);        //< jacobian J=det F

  // parameters
  static constexpr auto c1 = PARAM(0);   //< material parameter
  static constexpr auto c2 = PARAM(1);   //< material parameter
  static constexpr auto b = PARAM(2);   //< material parameter
  static constexpr auto d = PARAM(3);   //< material parameter

  static constexpr int nMaterialParameters = 4;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy density function, Psi_iso, in terms of the reduced invariants
  // function from chemo-electro-mechanical muscle model (Heidlauf 2013 "Modeling the chemo-electro-mechanical behaviour of... p.4)
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3)) + b/d * (pow(lambda, d) - INT(1)) - b*ln(lambda);
};

}  // namespace
}  // namespace
