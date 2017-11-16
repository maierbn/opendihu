#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t), R is from cellml file
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  OperatorSplitting::Godunov<
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::RegularFixed<1>,
        BasisFunction::Lagrange,
        Equation::Dynamic::Diffusion
      >
    >,
    TimeSteppingScheme::ExplicitEuler<
      CellmlAdapter
    >
  >
  operatorSplitting(settings);
       
  Computation computation(settings, operatorSplitting);
  computation.run();
  
  return EXIT_SUCCESS;
}