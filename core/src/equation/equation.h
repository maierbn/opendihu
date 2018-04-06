#pragma once

namespace Equation
{

/**
 * base class for equation
 */
class Equation
{
public:
  virtual ~Equation() {}
private:
};

/** No equation, can be used for testing I/O, because no solving is attempted
 */
class None : public Equation
{
public:
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = true;             ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = false;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = false;              ///< Equations of solid mechanics
};

}  // namespace
