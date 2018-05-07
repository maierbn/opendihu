#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 2D Mooney-Rivlin incompressible material, mixed formulation, Taylor-Hood elements, no static condensation
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::Mixed<
      BasisFunction::LagrangeOfOrder<1>,
      BasisFunction::LagrangeOfOrder<2>
    >,
    Quadrature::Mixed<
      Quadrature::Gauss<4>,  // low order (pressure)
      Quadrature::Gauss<3>   // high order (displacements)
    >,
    Equation::Static::MooneyRivlinIncompressible2D
  >
  problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
