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
  Control::PreciceAdapterVolumeCoupling<
    Control::Coupling<
      // multiple 1D fibers
      FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
        Control::MultipleInstances<                       // fibers
          OperatorSplitting::Strang<
            Control::MultipleInstances<
              TimeSteppingScheme::Heun<                   // fiber reaction term
                CellmlAdapter<
                  57, 71,  // nStates,nAlgebraics: 57,71 = Shorten, 4,9 = Hodgkin Huxley
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
      >,
      // 3D bidomain eq.
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
    >
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




