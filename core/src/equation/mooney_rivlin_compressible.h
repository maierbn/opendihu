#pragma once

#include "equation/static.h"

namespace Equation
{
namespace Static
{

/** compressible mooney rivlin u/p formulation, following Sussmann "A finite element formulation for nonlinear incompressible elastic and inelastic analysis"
 */
struct MooneyRivlinCompressible : public Static
{
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = false;            ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = false;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = true;               ///< Equations of solid mechanics
  static constexpr bool isIncompressible = false;              ///< Incompressible material
};

}  // namespace
}  // namespace
