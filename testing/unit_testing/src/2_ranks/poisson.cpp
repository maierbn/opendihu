#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"

TEST(PoissonTest, SerialEqualsParallelGlobalInput)
{
  // run serial problem
  std::string pythonConfig = R"(
# Poisson 2D,  3 x 4 = 12 nodes

nx = 2   # number of elements in x direction
ny = 3   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

# right hand side
rhs = [0.1, -0.1, 0.2, -0.2, 0.3, -0.3, 0.4, -0.4, 0.5, -0.5, 0.6, -0.6]

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": bc,
        "rightHandSide": rhs,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "relativeTolerance": 1e-15,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out22", "outputInterval": 1, "binary": False}
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
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Poisson
    >
  > ProblemType;

  ProblemType problemSerial(settings);
  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Poisson 2D,  3 x 4 = 12 nodes

nx = 2   # number of elements in x direction
ny = 3   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

# right hand side
rhs = [0.1, -0.1, 0.2, -0.2, 0.3, -0.3, 0.4, -0.4, 0.5, -0.5, 0.6, -0.6]

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": bc,
        "rightHandSide": rhs,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "relativeTolerance": 1e-15,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out22", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out22.py", "out22.0.py", "out22.1.py"};
  assertParallelEqualsSerialOutputFiles(outputFilesToCheck);

  nFails += ::testing::Test::HasFailure();
}

TEST(PoissonTest, SerialEqualsParallelLocalInput)
{
  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Poisson
    >
  > ProblemType;

  // run serial problem
  std::string pythonConfig = R"(
# Poisson 2D,  3 x 4 = 12 nodes

nx = 2   # number of elements in x direction
ny = 3   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

# right hand side
rhs = [0.1, -0.1, 0.2, -0.2, 0.3, -0.3, 0.4, -0.4, 0.5, -0.5, 0.6, -0.6]

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": bc,
        "rightHandSide": rhs,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "relativeTolerance": 1e-15,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out23", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  ProblemType problemSerial(settings);
  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Poisson 2D,  3 x 4 = 12 nodes (parallel)

# extract rankNo and nRanks from arguments
import sys
rankNo = (int)(sys.argv[-2])
nRanks = (int)(sys.argv[-1])
print("rankNo: {}, nRanks: {}".format(rankNo,nRanks))

nx = 2   # number of elements in x direction
ny = 3   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  if rankNo == 0:
    #bc[i] = np.sin(x*np.pi)
    bc[i] = i
  if rankNo == 1:
    i2 = (nx+1)*1 + i
    #bc[i2] = np.sin(x*np.pi)
    bc[i2] = 10*i

# right hand side
if rankNo == 0:
  rhs = [0.1, -0.1, 0.2,  -0.2, 0.3, -0.3]
  nElements = [nx, 2]
  physicalExtent = [2*nx, 4]

if rankNo == 1:
  rhs = [0.4, -0.4, 0.5,  -0.5, 0.6, -0.6]
  nElements = [nx, 1]
  physicalExtent = [2*nx, 2]

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": False,
        "dirichletBoundaryConditions": bc,
        "rightHandSide": rhs,
        "nElements": nElements,
        "nRanks": [1,2],
        "physicalExtent": physicalExtent,
        "relativeTolerance": 1e-15,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out23", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out23.py", "out23.0.py", "out23.1.py"};
  assertParallelEqualsSerialOutputFiles(outputFilesToCheck);

  nFails += ::testing::Test::HasFailure();
}
