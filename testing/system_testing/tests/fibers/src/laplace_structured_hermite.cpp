#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D Laplace equation
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  Postprocessing::StreamlineTracer<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::Hermite,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
      //Equation::None
    > 
  >
  problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
