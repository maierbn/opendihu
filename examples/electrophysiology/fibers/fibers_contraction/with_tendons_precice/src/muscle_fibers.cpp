#include <Python.h>
#include <iostream>
#include <cstdlib>

#include "easylogging++.h"

#include "opendihu.h"

//#define Shorten
#define HodgkinHuxley
//#define HodgkinHuxlexRazumova

// #define FiberDiffusionSolver TimeSteppingScheme::ImplicitEuler
// #define FiberDiffusionSolver TimeSteppingScheme::CrankNicolson
//template<typename T> using FiberDiffusionSolver = TimeSteppingScheme::ImplicitEuler<T>;
// template<typename T> using FiberDiffusionSolver = TimeSteppingScheme::CrankNicolson<T>;

////// ^^^^^^^^^^^^^^^^



#ifdef Shorten
#define N_STATES 57
#define N_ALGEBRAICS 1
#endif
#ifdef HodgkinHuxley
#define N_STATES 4
#define N_ALGEBRAICS 9
#endif
#ifdef HodgkinHuxlexRazumova
#define N_STATES 9
#define N_ALGEBRAICS 19
#endif

int main(int argc, char *argv[])
{
  // multiple fibers in arbitrary partitioning, coupled to dynamic nonlinear elasticity
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  
    Control::Coupling<
        FastMonodomainSolver<                               // a wrapper that improves performance of multidomain
            Control::MultipleInstances<                       // fibers
            OperatorSplitting::Strang<
                Control::MultipleInstances<
                TimeSteppingScheme::Heun<                   // fiber reaction term
                    CellmlAdapter<
                    N_STATES, N_ALGEBRAICS,                 // depends on the cellml model
                    FunctionSpace::FunctionSpace<
                        Mesh::StructuredDeformableOfDimension<1>,
                        BasisFunction::LagrangeOfOrder<1>
                    >
                    >
                >
                >,
                Control::MultipleInstances<
                TimeSteppingScheme::CrankNicolson<                       // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
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
        MuscleContractionSolver<>
        // OutputWriter::OutputSurface<
        // TimeSteppingScheme::StaticBidomainSolver<         // bidomain
        //     SpatialDiscretization::FiniteElementMethod<     // FEM for initial potential flow, fiber directions
        //         Mesh::StructuredDeformableOfDimension<3>,
        //         BasisFunction::LagrangeOfOrder<1>,
        //         Quadrature::Gauss<3>,
        //         Equation::Static::Laplace
        //     >,
        //     SpatialDiscretization::FiniteElementMethod<     // anisotropic diffusion
        //         Mesh::StructuredDeformableOfDimension<3>,
        //         BasisFunction::LagrangeOfOrder<1>,
        //         Quadrature::Gauss<5>,
        //         Equation::Dynamic::DirectionalDiffusion
        //     >
        // >
        // >
    > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




