#pragma once

#include "equation/dynamic.h"

namespace Equation
{
namespace Dynamic
{

/** u_t - Δu = 0
  */
struct IsotropicDiffusion : public Dynamic
{
  static constexpr bool usesTimeStepping = true;               ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = true;             ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = true;                         ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = false;              ///< Equations of solid mechanics
};

/** u_t - ∇•A∇u = 0
  */
struct AnisotropicDiffusion : public Dynamic
{
  static constexpr bool usesTimeStepping = true;               ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = false;            ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = true;  ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = true;                         ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = false;              ///< Equations of solid mechanics
};

/** u_t - ∇•A(v)∇u = 0 where the diffusion tensor, A, depends upon a direction field, v
  */
struct DirectionalDiffusion : public AnisotropicDiffusion
{
};

}  // namespace
}  // namespace
