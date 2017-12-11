#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D diffusion equation du/dt = c du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::RegularFixed<2>,
      BasisFunction::Lagrange,
      Equation::Dynamic::Diffusion
    >
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  return EXIT_SUCCESS;
}