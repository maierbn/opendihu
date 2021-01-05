#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D Laplace equation 0 = du^2/dx^2 + du^2/dy^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();
  
  return EXIT_SUCCESS;
}
