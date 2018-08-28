#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"

TEST(LaplaceTest, Test1)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[-1] = 0.0

config = {
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "nElements": n,
    "physicalExtent": 4.0,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > problem(settings);
  
  problem.run();
}
