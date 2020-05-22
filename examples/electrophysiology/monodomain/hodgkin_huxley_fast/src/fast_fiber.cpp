#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // A single Hodgkin-Huxley fiber using the FastMonodomainSolver.
  // The purpose of this file is to compare with not_fast_fiber.cpp
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem with FastMonodomainSolver
  TimeSteppingScheme::RepeatedCall<
    FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
      Control::MultipleInstances<                       // fibers
        OperatorSplitting::Strang<
          Control::MultipleInstances<
            TimeSteppingScheme::Heun<                   // fiber reaction term
              CellmlAdapter<
                4, 9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
                FunctionSpace::FunctionSpace<
                  Mesh::StructuredDeformableOfDimension<1>,
                  BasisFunction::LagrangeOfOrder<1>
                >
              >
            >
          >,
          Control::MultipleInstances<
            TimeSteppingScheme::ImplicitEuler<          // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
              SpatialDiscretization::FiniteElementMethod<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>,
                Quadrature::Gauss<2>,
                Equation::Dynamic::IsotropicDiffusion
              >
            >
          >
        >
      >
    >
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




