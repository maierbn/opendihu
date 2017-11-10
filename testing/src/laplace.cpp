#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"

TEST(LaplaceTest, Matrix)
{
  // initialize everything, handle arguments and parse settings from input file

  std::string pythonConfig = R"(
# Laplace 1D

n = 40


# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtend": 4.0,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  },
  "OutputWriter" : [
    {"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
    {"format": "Python", "filename": "p"}
  ]
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::RegularFixed<1>,
    BasisFunction::Lagrange,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
}

