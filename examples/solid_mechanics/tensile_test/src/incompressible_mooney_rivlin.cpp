#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

// define material
struct Material : Equation::SolidMechanics::HyperelasticityBase
{
  static constexpr bool isIncompressible = true;       //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = false;    //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress
  
  // material parameters
  static constexpr auto c1 = PARAM(0);     //< material parameter
  static constexpr auto c2 = PARAM(1);     //< material parameter

  static constexpr int nMaterialParameters = 2;  //< number of material parameters

  //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3));

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

  //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C.
  //! it must only depend on variables C11, C12, C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);
};

int main(int argc, char *argv[])
{
  // Solves nonlinear hyperelasticity (Mooney-Rivlin) using the built-in solver
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  
  // define problem
  SpatialDiscretization::HyperelasticitySolver<Material> problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




