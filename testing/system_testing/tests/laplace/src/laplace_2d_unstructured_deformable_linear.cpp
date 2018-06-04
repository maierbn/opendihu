#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 2D laplace equation du^2/dx^2 = 0
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  

  SpatialDiscretization::FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
