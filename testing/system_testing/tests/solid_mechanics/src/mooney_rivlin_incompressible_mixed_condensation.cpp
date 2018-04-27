#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D Mooney-Rivlin incompressible material, mixed formulation with constant pressure and static condensation
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::Mixed<
      BasisFunction::CompletePolynomialOfDimensionAndOrder<3,0>,
      BasisFunction::LagrangeOfOrder<2>
    >,
    Quadrature::Mixed<
      Quadrature::Gauss<2>,
      Quadrature::Gauss<3>
    >,
    Equation::Static::MooneyRivlinIncompressible3D
  >
  problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
