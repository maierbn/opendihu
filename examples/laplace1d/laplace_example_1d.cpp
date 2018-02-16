#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D Laplace equation 0 = du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange<1>,
    Integrator::Gauss<2>,
    Equation::Static::Laplace
  > problem(settings);
  
  problem.run();
  
  /*Computation computation(settings, equationDiscretized);
  computation.run();*/
  
  return EXIT_SUCCESS;
}