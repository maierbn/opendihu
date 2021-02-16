#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  Control::Coupling<
    // prescribed activation
    PrescribedValues<
      FunctionSpace::FunctionSpace<
        Mesh::CompositeOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>
      >
    >,
    // quasi-static mechanics solver
    MuscleContractionSolver<
      Mesh::CompositeOfDimension<3>
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
