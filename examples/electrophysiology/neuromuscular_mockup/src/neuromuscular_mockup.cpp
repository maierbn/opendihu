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
  
  // define helper function space for various activation signals, this is actually a vector space
  using HelperFunctionSpace = FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>>;
  
  // define overall structure of solvers
  Control::Coupling<
    // motoneuron input solver
    Control::Coupling<  
      // mapping muscle spindles -> motor neuron signals
      Control::MapDofs<
        HelperFunctionSpace,
        // muscle spindles solver
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,
            HelperFunctionSpace
          >
        >
      >,
      // processed Golgi tendon sensor solver
      Control::Coupling<
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
        >
      >
    >,
    // motoneuron with electro-mechanics
    Control::Coupling<
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
        
        // map from λ in the 3D mesh to golgi tendon organs
        Control::MapDofs<
          HelperFunctionSpace,
          
          // map from motoneuronMesh to stimulated nodes
          Control::MapDofs<
            HelperFunctionSpace,
            
            // electro-mechanics solver
            Control::Coupling
            <
              // electrophysiology solver (mockup)
              Control::MultipleInstances<
                PrescribedValues<
                  FunctionSpace::FunctionSpace<
                    Mesh::StructuredDeformableOfDimension<3>,
                    BasisFunction::LagrangeOfOrder<1>
                  >
                >
              >,
              // mechanics solver
              MuscleContractionSolver<
                Mesh::CompositeOfDimension<3>
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




