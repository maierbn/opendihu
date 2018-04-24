#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 2D diffusion equation du/dt = c du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::Hermite,
      Quadrature::Gauss<3>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
