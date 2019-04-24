#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D Laplace equation 0 = du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > problem(settings);
  
  problem.run();
  
  /* 
  equationDiscretized.run();*/
  
  return EXIT_SUCCESS;
}
