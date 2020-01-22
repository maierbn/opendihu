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
    TimeSteppingScheme::ImplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<1>,
        BasisFunction::LagrangeOfOrder<>,
        Quadrature::None,
        Equation::Dynamic::IsotropicDiffusion
        >
    >
  > problem(settings);

  // run the simulation
  problem.run();

  return EXIT_SUCCESS;
}
