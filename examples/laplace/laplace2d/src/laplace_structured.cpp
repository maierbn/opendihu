#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 2D Laplace equation 0 = du^2/dx^2 + du^2/dy^2
  
  // initialize and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define the tree of solvers (here only one FEM solver)
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  // run the simulation
  equationDiscretized.run();
  
  return EXIT_SUCCESS;
}
