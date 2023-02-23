#include <iostream>
#include <cstdlib>

#include "opendihu.h"


typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> HelperFunctionSpace;

int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  TimeSteppingScheme::Heun<
    CellmlAdapter<
      8,23,
      HelperFunctionSpace
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
