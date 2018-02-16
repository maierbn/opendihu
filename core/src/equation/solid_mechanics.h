#pragma once

#include "equation/static.h"

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace Equation
{
namespace Static
{

/** Hyperelastic solid mechanics formulation with arbitrary strain energy density function
 */
class SolidMechanics : public Static
{
public:
     
  // define helper variables for strain energy density function
  // invariants
  static constexpr auto I1 = VAR(0);
  static constexpr auto I2 = VAR(1);
  static constexpr auto I3 = VAR(2);
    
  // parameters
  static constexpr auto c0 = PARAM(0);
  static constexpr auto c1 = PARAM(1);

  static const auto constexpr strainEnergyDensityFunction 
    = c0*(I1 - INT(3)) + c1*(I2 - INT(3)) + sin(I3) + pow(I3,INT(5));
  
private:
};

}  // namespace
}  // namespace
