#include <iostream>
#include <cstdlib>

#include "opendihu.h"

// define helper function space for various activation signals, this is actually a vector space
typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> HelperFunctionSpace;
  
int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  Control::Coupling<
    Control::MultipleCoupling<
      // muscle spindles solver
      Control::MapDofs<
        HelperFunctionSpace,
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            8,23,
            HelperFunctionSpace
          >
        >
      >,
      // golgi tendon organ
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          4,9,
          HelperFunctionSpace
        >
      >,
      // motor neurons
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          6,14,
          HelperFunctionSpace
        >
      >,
      // Fast monodomain: 0D / 1D + 3D bidomain
      Control::MapDofs<
        HelperFunctionSpace,
        Control::MapDofs<
          HelperFunctionSpace,
          Control::Coupling<
            // Fast monodomain: 0D / 1D
            FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
              Control::MultipleInstances<                       // fibers
                OperatorSplitting::Strang<
                  Control::MultipleInstances<
                    TimeSteppingScheme::Heun<                   // fiber reaction term
                      CellmlAdapter<
                        9, 19,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley. 9.19 = HH-R
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
            // 3D bidomain
            OutputWriter::OutputSurface<
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
          >
        >
      >,
      // Fast monodomain: 0D / 1D + 3D bidomain
      // Control::MapDofs<
      //   HelperFunctionSpace,
        Control::Coupling<
          // Fast monodomain: 0D / 1D
          FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
            Control::MultipleInstances<                       // fibers
              OperatorSplitting::Strang<
                Control::MultipleInstances<
                  TimeSteppingScheme::Heun<                   // fiber reaction term
                    CellmlAdapter<
                      9, 19,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley. 9.19 = HH-R
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
          // 3D bidomain
          OutputWriter::OutputSurface<
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
        >
      // >
    >,
    // 2x mechanics TODO make sure that no data is copied
    Control::Coupling<
      MuscleContractionSolver<
        Mesh::StructuredDeformableOfDimension<3>
      >,
      MuscleContractionSolver<
        Mesh::StructuredDeformableOfDimension<3>
      >
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
