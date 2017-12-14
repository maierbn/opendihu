#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "equation/diffusion.h"

namespace SpatialDiscretization
{
  
TEST(DiffusionTest, Compiles)
{
  std::string pythonConfig = R"(
# Diffusion 1D
n = 5
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 1.0,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtend": 4.0,
      "relativeTolerance": 1e-15,
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "Python", "filename": "p", "outputInterval": 1}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  TimeSteppingScheme::ExplicitEuler<
    FiniteElementMethod<
      Mesh::RegularFixed<1>,
      BasisFunction::Lagrange<>,
      Integrator::None,
      Equation::Dynamic::Diffusion
    >
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
}

} // namespace
