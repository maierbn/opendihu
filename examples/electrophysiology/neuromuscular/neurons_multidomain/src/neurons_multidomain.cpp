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
    // prescribed activation
    PrescribedValues<
      FunctionSpace::FunctionSpace<
        Mesh::CompositeOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>
      >
    >,
    // quasi-static mechanics solver for "precontraction"
    MuscleContractionSolver<
      Mesh::CompositeOfDimension<3>
    >,
    // static trans-iso material for "prestretch"
    SpatialDiscretization::HyperelasticitySolver<
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D, true, Mesh::CompositeOfDimension<3>
    >,
    // actual simlation
    
    Control::MultipleCoupling<
      // mapping muscle spindles -> motor neuron signals
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
      // map from λ in the 3D mesh to muscle spindles input (here: just prescribe sensor values)
      Control::MapDofs<
        HelperFunctionSpace,
        
        // map from λ in the 3D mesh to golgi tendon organs (here: just prescribe sensor values)
        Control::MapDofs<
          HelperFunctionSpace,
          
          // map from motoneuronMesh to stimulated nodes
          Control::MapDofs<
            HelperFunctionSpace,
            
            Control::Coupling<
              // Multidomain, Strang splitting of CellmlAdapter and MultidomainWithFatSolver
              OperatorSplitting::Strang<
                Control::MultipleInstances<
                  TimeSteppingScheme::Heun<
                    CellmlAdapter<
                      9,19,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley
                      FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
                    >  
                  >
                >,
                OutputWriter::OutputSurface<
                  TimeSteppingScheme::MultidomainWithFatSolver<       // multidomain
                    SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
                      MeshType,
                      BasisFunction::LagrangeOfOrder<1>,
                      Quadrature::Gauss<3>,
                      Equation::Static::Laplace
                    >,
                    SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
                      MeshType,
                      BasisFunction::LagrangeOfOrder<1>,
                      Quadrature::Gauss<3>,
                      Equation::Dynamic::DirectionalDiffusion
                    >,
                    SpatialDiscretization::FiniteElementMethod<       // isotropic diffusion in fat layer
                      MeshType,
                      BasisFunction::LagrangeOfOrder<1>,
                      Quadrature::Gauss<3>,
                      Equation::Dynamic::IsotropicDiffusion
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
