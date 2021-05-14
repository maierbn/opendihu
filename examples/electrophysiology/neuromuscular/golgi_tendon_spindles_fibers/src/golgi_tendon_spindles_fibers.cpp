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

  Control::MultipleCoupling<
    // Term1: prescribed activation
    PrescribedValues<
      FunctionSpace::FunctionSpace<
        Mesh::CompositeOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>
      >
    >,
    // Term2: quasi-static mechanics solver for "precontraction"
    MuscleContractionSolver<
      Mesh::CompositeOfDimension<3>
    >,
    // Term3: static trans-iso material for "prestretch"
    SpatialDiscretization::HyperelasticitySolver<
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D, true, Mesh::CompositeOfDimension<3>
    >,

    // Term4: actual simlation
    Control::MultipleCoupling<
      // mapping muscle spindles output -> motor neuron input
      Control::MapDofs<
        HelperFunctionSpace,
        // muscle spindles solver
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            8,23,
            HelperFunctionSpace
          >
        >
      >,
      // mapping Golgi tendon organs -> interneurons
      Control::MapDofs<
        HelperFunctionSpace,
        // Golgi tendon organs solver
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,
            HelperFunctionSpace
          >
        >
      >,
      // mapping interneurons -> input for motor neurons
      Control::MapDofs<
        HelperFunctionSpace,
        // interneurons solver
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,
            HelperFunctionSpace
          >
        >
      >,
      // mapping motor neuron signals + cortical input to actual inputs
      Control::MapDofs<
        HelperFunctionSpace,
        // motoneuron solver
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            6,14,  // nStates,nAlgebraics
            HelperFunctionSpace
          >
        >
      >,
      // map from λ in the 3D mesh to muscle spindles input
      Control::MapDofs<
        HelperFunctionSpace,
      
        // map from λ in the 3D mesh to golgi tendon organs (here: just prescribe sensor values)
        Control::MapDofs<
          HelperFunctionSpace,
            
          // map from motoneuronMesh to stimulated nodes
          Control::MapDofs<
            HelperFunctionSpace,
                    
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
              >,
              
              // muscle contraction
              MuscleContractionSolver<
                Mesh::CompositeOfDimension<3>
              >
            >
          >
        >
      >
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
