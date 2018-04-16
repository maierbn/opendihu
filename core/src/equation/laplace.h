#pragma once

#include "equation/static.h"

namespace Equation
{
namespace Static
{

/** Δu = 0
  */
struct Laplace : public Static
{
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = true;             ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = false;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = false;              ///< Equations of solid mechanics
};

/** ∇•A∇u = 0
  */
struct GeneralizedLaplace : public Static
{
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = false;            ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = true;  ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = false;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = false;              ///< Equations of solid mechanics
};

}  // namespace
}  // namespace
