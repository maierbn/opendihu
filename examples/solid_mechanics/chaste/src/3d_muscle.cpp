#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // adapter to chaste, solves nonlinear hyperelasticity (Mooney-Rivlin) using chaste with the muscle geometry
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  TimeSteppingScheme::QuasiStaticNonlinearElasticitySolverChaste<3> problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




