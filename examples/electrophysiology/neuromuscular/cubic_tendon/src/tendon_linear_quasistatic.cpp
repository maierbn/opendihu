#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  //SpatialDiscretization::HyperelasticitySolver<
  TimeSteppingScheme::QuasistaticHyperelasticitySolver<
    Equation::SolidMechanics::SaintVenantKirchhoff
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




