#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // Solves nonlinear hyperelasticity (Mooney-Rivlin) using the built-in solver with the muscle geometry
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  TimeSteppingScheme::QuasistaticHyperelasticitySolver<
        Equation::SolidMechanics::HyperelasticTendon
        //Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D
  
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




