#pragma once

#include "equation/static.h"

namespace Equation
{
namespace Static
{

/** compressible mooney rivlin u/p formulation, following Sussmann "A finite element formulation for nonlinear incompressible elastic and inelastic analysis"
 */
struct CompressibleMooneyRivlin : public Static
{
  static const bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static const bool hasLaplaceOperator = false;            ///< Equations that include Δu
  static const bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static const bool hasRhs = false;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static const bool isSolidMechanics = true;               ///< Equations of solid mechanics

};

}  // namespace
}  // namespace
