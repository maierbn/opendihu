#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D multidomain coupled with contraction
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;

  Control::Coupling<
    // prestretch, static trans-iso
    SpatialDiscretization::HyperelasticitySolver<
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D, true, Mesh::CompositeOfDimension<3>
    >,
    Dummy
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
