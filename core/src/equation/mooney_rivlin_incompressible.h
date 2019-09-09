#pragma once

#include "equation/static.h"

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace Equation
{
namespace SolidMechanics
{

/** Hyperelastic solid mechanics formulation with arbitrary strain energy density function
 */
struct MooneyRivlinIncompressible3D
{
  // define helper variables for strain energy density function
  // reduced invariants
  static constexpr auto Ibar1 = VAR(0);    //< 1st reduced invariant
  static constexpr auto Ibar2 = VAR(1);    //< 2nd reduced invariant
  static constexpr auto J = VAR(2);        //< jacobian J=det F

  // parameters
  static constexpr auto c1 = PARAM(0);   //< material parameter
  static constexpr auto c2 = PARAM(1);   //< material parameter

  static constexpr int nMaterialParameters = 2;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy density function, Psi_iso, in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3));
};

}  // namespace
}  // namespace
