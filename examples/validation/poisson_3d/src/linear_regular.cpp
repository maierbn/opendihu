#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D Poisson equation f = du^2/dx^2
    
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Poisson
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
