#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "../utility.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

TEST(SolidMechanicsTest, Test0)
{/*
  std::string pythonConfig = R"(
# solid mechanics (3D)

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": [2,2,2],   # 8 elements
    "physicalExtent": [4.0,4.0,4.0],
    "dirichletBoundaryConditions": {0:1.0},
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "PythonFile", "filename" : "out_txt", "binary" : False},
    ]
  }
}
)";
  DihuContext settings(argc, argv, pythonConfig);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::Mixed<
      BasisFunction::LagrangeOfOrder<1>,
      BasisFunction::LagrangeOfOrder<2>
    >,
    Quadrature::Mixed<
      Quadrature::Gauss<2>,
      Quadrature::Gauss<3>
    >,
    Equation::Static::SolidMechanics
  > problem(settings);
  
  problem.run();*/
}


