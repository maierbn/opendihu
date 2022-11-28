#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

// define custom material
struct Material : Equation::SolidMechanics::HyperelasticityBase
{
  static constexpr bool isIncompressible = false;      //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection = false;     //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar5, this means it is an anisotropic material
  static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress, this is only possible with MuscleContractionSolver
  
  // material parameters that will be set by the python settings, 
  // these constants can be used in all the strain energy functions below
  static constexpr auto c1 = PARAM(0);     //< material parameter
  static constexpr auto c2 = PARAM(1);     //< material parameter
  static constexpr auto c = PARAM(2);      //< material parameter

  static constexpr int nMaterialParameters = 3;  //< number of material parameters, this has to match the number of defined material parameters

  //! The isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants.
  //! Normally, this is used for incompressible materials, then set isIncompressible to true.
  //! Set usesFiberDirection to true if you use lambda, Ibar4 or Ibar5.
  //! It should only use the variables Ibar1, Ibar2, Ibar4, Ibar5, lambda (lambda is an alias for sqrt(Ibar4))
  static const auto constexpr strainEnergyDensityFunctionIsochoric = INT(0);

  //! The volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false), J = det(F)
  //! It should only use the variable J.
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! Coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions, I1 = tr(C), I2 = 1/2 (tr(C)^2 - tr(C^2)), I3 = det(C) = J^2
  //! It should use the variables I1, I2, I3.
  static constexpr auto d = INT(2)*(c1 + INT(2)*c2);
  static const auto constexpr strainEnergyDensityFunctionCoupled = c*pow(sqrt(I3) - INT(1), INT(2)) - d*ln(sqrt(I3)) + c1*(I1 - INT(3)) + c2*(I2 - INT(3));

  //! Another coupled form of the strain energy function, Ψ(C), dependend on right Cauchy Green tensor C=(C11,...,C33)  and fiber direction a=(a1,a2,a3).
  //! Set usesFiberDirection to true if you use this.
  //! It should only use the variables C11, C12, C13, C22, C23, C33, a1, a2, a3, I4 (I4 is an alias for a0•C a0) (I5 can be computed manually from C and a).
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);
};

int main(int argc, char *argv[])
{
  // Solves nonlinear hyperelasticity (Mooney-Rivlin) using the built-in solver with the muscle geometry
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  SpatialDiscretization::HyperelasticitySolver<
    Material
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




