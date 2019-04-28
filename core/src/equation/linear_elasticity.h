#pragma once

#include "equation/static.h"

namespace Equation
{
namespace Static
{

/** ∫σ:ε dx = rhs
  */
struct LinearElasticity : public Static
{
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = true;             ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = true;                         ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = false;              ///< Equations of solid mechanics
};

}  // namespace
}  // namespace
