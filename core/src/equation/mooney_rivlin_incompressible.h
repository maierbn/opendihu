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
struct MooneyRivlinIncompressible3D : public Static
{
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = false;            ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = true;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = true;               ///< Equations of solid mechanics
  static constexpr bool isIncompressible = true;              ///< Incompressible material
     
  // define helper variables for strain energy density function
  // reduced invariants
  static constexpr auto Ibar1 = VAR(0);    //< 1st reduced invariant
  static constexpr auto Ibar2 = VAR(1);    //< 2nd reduced invariant
  static constexpr auto J = VAR(2);     //< jacobian J=det F
    
  // parameters
  static constexpr auto c0 = PARAM(0);   //< material parameter
  static constexpr auto c1 = PARAM(1);   //< material parameter
  static constexpr auto kappa = PARAM(2);   //< penalty parameter, kappa->infinity is incompressible
  
  static constexpr int nMaterialParameters = 3;  //< number of material parameters
  
  
  //! the isochoric part of the decoupled strain energy density function, Psi_iso, in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric 
    = c0*(Ibar1 - INT(3)) + c1*(Ibar2 - INT(3));
    
  //! the volumetric part of the strain energy density function, numerical constant (kappa) * penalty function
  static const auto constexpr strainEnergyDensityFunctionVolumetric
    = kappa/INT(2) * pow(J - INT(1), INT(2));
};

/** Hyperelastic solid mechanics formulation with arbitrary strain energy density function
 */
struct MooneyRivlinIncompressible2D : public Static
{
  static constexpr bool usesTimeStepping = false;              ///< Equation of the form L = u_t
  static constexpr bool hasLaplaceOperator = false;            ///< Equations that include Δu
  static constexpr bool hasGeneralizedLaplaceOperator = false; ///< Equations that include ∇•(A∇u)
  static constexpr bool hasRhs = true;                        ///< Equations that can have a non-zero rhs (Lu = f)
  static constexpr bool isSolidMechanics = true;               ///< Equations of solid mechanics
  static constexpr bool isIncompressible = true;              ///< Incompressible material
     
  // define helper variables for strain energy density function
  // reduced invariants
  static constexpr auto Ibar1 = VAR(0);    //< 1st reduced invariant
  static constexpr auto Ibar2 = VAR(1);    //< 2nd reduced invariant
  static constexpr auto J = VAR(2);     //< jacobian J=det F
    
  // parameters
  static constexpr auto c0 = PARAM(0);   //< material parameter
  static constexpr auto c1 = PARAM(1);   //< material parameter
  static constexpr auto kappa = PARAM(2);   //< penalty parameter, kappa->infinity is incompressible
  
  static constexpr int nMaterialParameters = 3;  //< number of material parameters
  
  
  //! the isochoric part of the decoupled strain energy density function, Psi_iso, in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric 
    = c0*(Ibar1 - INT(2)) + c1*(Ibar2 - INT(2));
    
  //! the volumetric part of the strain energy density function, numerical constant (kappa) * penalty function
  static const auto constexpr strainEnergyDensityFunctionVolumetric
    = kappa/INT(2) * pow(J - INT(1), INT(2));
};
}  // namespace
}  // namespace
