#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 2D multidomain: implicit Euler, FEM
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  OperatorSplitting::Godunov<
    Control::MultipleInstances<
      TimeSteppingScheme::ExplicitEuler<     // Hodgkin-Huxley
        CellmlAdapter<4>
      >,
    >,
    TimeSteppingScheme::MultidomainSolver<     // multidomain
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
