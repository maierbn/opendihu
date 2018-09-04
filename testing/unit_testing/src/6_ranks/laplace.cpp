#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"

TEST(LaplaceTest, SerialEqualsParallelStructured1DLinear)
{
  // run serial problem
  std::string pythonConfig = R"(
import numpy as np

nx = 2   # number of elements in x direction
ny = 3   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(2*nx+1)):
  x = i/(2*nx+1.)
  bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (2*nx+1)*ny + i
  bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "DirichletBoundaryCondition": bc,
        "relativeTolerance": 1e-15,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<2>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();

  LOG(INFO) << " =================== run parallel problem ================= ";

  // run parallel problem
  std::string pythonConfig2 = R"(
import numpy as np

nx = 2   # number of elements in x direction
ny = 3   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(2*nx+1)):
  x = i/(2*nx+1.)
  bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (2*nx+1)*ny + i
  bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "DirichletBoundaryCondition": bc,
        "relativeTolerance": 1e-15,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
        ]
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
