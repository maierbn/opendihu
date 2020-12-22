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
  Control::Coupling<
    Control::MultipleInstances<
      PrescribedValues<               // fibers with prescribed values
        FunctionSpace::FunctionSpace<
          Mesh::StructuredDeformableOfDimension<1>,
          BasisFunction::LagrangeOfOrder<1>
        >
      >
    >,
    TimeSteppingScheme::StaticBidomainSolver<           // Bidomain solver
      SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<3>,
        Equation::Static::Laplace
      >,
      SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<5>,
        Equation::Dynamic::DirectionalDiffusion
      >
    >
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




