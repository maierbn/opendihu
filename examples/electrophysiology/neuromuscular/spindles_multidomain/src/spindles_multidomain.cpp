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
      // motoneuron solver
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          6,14,  // nStates,nAlgebraics
          HelperFunctionSpace
        >
      >,
      // map from λ in the 3D mesh to muscle spindles input
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
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
