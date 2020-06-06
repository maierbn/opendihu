#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // Solves linear solid mechanics using the built-in solver
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  
  // define problem
  SpatialDiscretization::FiniteElementMethod<       // linear elasticity
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::LagrangeOfOrder<2>,        // quadratic mesh, the same as for nonlinear mechanics, such that BCs can be reused
      Quadrature::Gauss<3>,
      Equation::Static::LinearElasticity
    > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




