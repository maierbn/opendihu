#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D Poisson equation f = du^2/dx^2
    
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Equation::Static::Poisson
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  return EXIT_SUCCESS;
}
