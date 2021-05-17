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
  
  Control::MultipleCoupling<
    // prescribed activation
    PrescribedValues<
      FunctionSpace::FunctionSpace<
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>
      >
    >,
    // quasi-static mechanics solver for "precontraction"
    MuscleContractionSolver<
      Mesh::StructuredDeformableOfDimension<3>
    >,
    // static trans-iso material for "prestretch"
    SpatialDiscretization::HyperelasticitySolver<
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D, true, Mesh::StructuredDeformableOfDimension<3>
    >,
    // actual simlation
    Control::PreciceAdapterVolumeCoupling<
      MuscleContractionSolver<>
    >
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




