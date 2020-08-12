#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D contraction
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
    
  Control::Coupling
  <
    // electrophysiology solver (mockup)
    Control::MultipleInstances<
      PrescribedValues<
        FunctionSpace::FunctionSpace<
          Mesh::StructuredDeformableOfDimension<3>,
          BasisFunction::LagrangeOfOrder<1>
        >
      >
    >,
    // mechanics solver
    MuscleContractionSolver<
      Mesh::CompositeOfDimension<3>
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
