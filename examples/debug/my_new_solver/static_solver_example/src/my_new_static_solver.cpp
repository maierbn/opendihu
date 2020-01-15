#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 2D Laplace quation with linear ansatz functions, cΔu = 0
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define the problem type, the nested solvers
  MyNewStaticSolver<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
    >
  > problem(settings);
    
  // run the simulation
  problem.run();
    
  return EXIT_SUCCESS;
}
