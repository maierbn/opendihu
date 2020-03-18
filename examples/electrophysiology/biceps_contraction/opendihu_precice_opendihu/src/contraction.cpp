#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // multiple fibers in arbitrary partitioning, coupled to dynamic nonlinear elasticity
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  PreciceAdapter::MuscleContraction<
    MuscleContractionSolver
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




