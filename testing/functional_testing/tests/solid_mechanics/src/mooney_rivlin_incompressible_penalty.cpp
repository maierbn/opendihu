#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D Mooney-Rivlin incompressible material, penalty formulation
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::MooneyRivlinIncompressible
  >
  problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
