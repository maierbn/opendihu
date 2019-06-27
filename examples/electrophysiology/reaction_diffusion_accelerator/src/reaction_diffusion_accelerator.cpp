#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // multiple fibers in arbitrary partitioning, coupled to bidomain equation for EMG
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  Control::Coupling<
    Control::ReactionDiffusionAccelerator< 57, 1>,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
    TimeSteppingScheme::StaticBidomainSolver<              // bidomain
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




