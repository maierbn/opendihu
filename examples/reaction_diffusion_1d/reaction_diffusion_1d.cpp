#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t)
  
  // initialize everything, handle arguments and parse settings from input file
  DiHuContext settings(argc, argv);

  OperatorSplitting::Godunov<
    TimeSteppingScheme::ExplicitEuler<
      Equation::Dynamic::ReactionDiffusion::ReactionTerm
    >,
    TimeSteppingScheme::CrankNicholson<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::Regular<1>,
        BasisFunction::Lagrange<1>,
        Equation::Dynamic::ReactionDiffusion::DiffusionTerm
      >
    >
  >
  operatorSplitting(settings);
                                                        
  Computation computation(settings, operatorSplitting);
  computation.run();
                                                      
  
  return EXIT_SUCCESS;
}