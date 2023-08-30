#include <iostream>
#include <cstdlib>

#include "opendihu.h"

typedef Mesh::StructuredDeformableOfDimension<3> MeshType;

// define helper function space for various activation signals, this is actually a vector space
typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> HelperFunctionSpace;
  
int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  // electrophysiology + mechanics
  Control::Coupling<
    
    // fibers emg, i.e., electrophysiology with multiple fibers, monodomain and bidomain equations
    Control::Coupling<
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
      OutputWriter::OutputSurface<
        TimeSteppingScheme::StaticBidomainSolver<              // bidomain
          SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
            Mesh::CompositeOfDimension<3>,
            BasisFunction::LagrangeOfOrder<1>,
            Quadrature::Gauss<3>,
            Equation::Static::Laplace
          >,
          SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
            Mesh::CompositeOfDimension<3>,
            BasisFunction::LagrangeOfOrder<1>,
            Quadrature::Gauss<5>,
            Equation::Dynamic::DirectionalDiffusion
          >
        >
      >
    >,
    
    // muscle contraction
    MuscleContractionSolver<
      Mesh::CompositeOfDimension<3>
    >
  > problem(settings);

  problem.run();
  
  return EXIT_SUCCESS;
}
