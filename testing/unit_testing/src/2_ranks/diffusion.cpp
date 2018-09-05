#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"


// 2D structured deformable
TEST(DiffusionTest, SerialEqualsParallelDeformable2DLinear)
{
  // run serial problem
  std::string pythonConfig = R"(
# Diffusion 2D

nx = 5   # number of elements in x direction
ny = 8   # number of elements in y direction

# initial values
iv = {}
iv[22] = 5.
iv[23] = 4.
iv[27] = 4.
iv[28] = 3.

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "ExplicitEuler": {
        "initialValues": iv,
        "numberTimeSteps": 5,
        "endTime": 0.1,
        "FiniteElementMethod": {
          "inputMeshIsGlobal": True,
          "nElements": [nx, ny],
          "physicalExtent": [2*nx, 2*ny],
          "relativeTolerance": 1e-15,
          "OutputWriter" : [
            {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
          ]
        }
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef Control::MultipleInstances<
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<2>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<2>,
        Equation::Static::Laplace
      >
    >
  > ProblemType;

  ProblemType problemSerial(settings);
  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Diffusion 2D

nx = 5   # number of elements in x direction
ny = 8   # number of elements in y direction

# initial values
iv = {}
iv[22] = 5.
iv[23] = 4.
iv[27] = 4.
iv[28] = 3.

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "ExplicitEuler": {
        "initialValues": iv,
        "numberTimeSteps": 5,
        "endTime": 0.1,
        "FiniteElementMethod": {
          "inputMeshIsGlobal": True,
          "nElements": [nx, ny],
          "physicalExtent": [2*nx, 2*ny],
          "relativeTolerance": 1e-15,
          "OutputWriter" : [
            {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
          ]
        }
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out.py", "out.0.py", "out.1.py"};
  assertParallelEqualsSerialOutputFiles(outputFilesToCheck);

  nFails += ::testing::Test::HasFailure();
}
