#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  Control::MultipleCoupling<
    // prescribed activation
    PrescribedValues<
      FunctionSpace::FunctionSpace<
        Mesh::CompositeOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>
      >
    >,
    // quasi-static mechanics solver for "precontraction"
    MuscleContractionSolver<
      Mesh::CompositeOfDimension<3>
    >,
    // static trans-iso material for "prestretch"
    SpatialDiscretization::HyperelasticitySolver<
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D, true, Mesh::CompositeOfDimension<3>
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
