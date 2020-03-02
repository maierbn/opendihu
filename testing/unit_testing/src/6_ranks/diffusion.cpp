#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"

// 2D regular fixed
TEST(DiffusionTest, SerialEqualsParallelGlobal)
{
  typedef Control::MultipleInstances<
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<2>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<2>,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  > ProblemType;

  // run global settings problem
  std::string pythonConfig = R"(
# Diffusion 2D

import sys
rankNo = (int)(sys.argv[-2])
nRanks = (int)(sys.argv[-1])

print(rankNo,nRanks)

nx = 13   # number of elements in x direction
ny = 8   # number of elements in y direction

# 14*9=126 nodes

# initial values
iv = {}
iv[62] = 5.
iv[63] = 4.5
iv[64] = 4.
iv[79] = 4.5
iv[80] = 4.
iv[81] = 3.5

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1,2,3,4,5],
      "ExplicitEuler": {
        "initialValues": iv,
        "numberTimeSteps": 5,
        "endTime": 0.1,
        "FiniteElementMethod": {
          "inputMeshIsGlobal": True,
          "nElements": [nx, ny],
          "physicalExtent": [2*nx, 2*ny],
          "relativeTolerance": 1e-15,
          "dirichletBoundaryConditions": {},
          "neumannBoundaryConditions": [],
        },
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out0", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  int ownRankNo = settings.ownRankNo();

  ProblemType problemGlobalSettings(settings);
  problemGlobalSettings.run();

  // run problem with global settings with only 1 rank
  std::string pythonConfig0 = R"(
# Diffusion 2D

import sys
rankNo = (int)(sys.argv[-2])
nRanks = (int)(sys.argv[-1])

print(rankNo,nRanks)

nx = 13   # number of elements in x direction
ny = 8   # number of elements in y direction

# 14*9=126 nodes

# initial values
iv = {}
iv[62] = 5.
iv[63] = 4.5
iv[64] = 4.
iv[79] = 4.5
iv[80] = 4.
iv[81] = 3.5

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [5],
      "ExplicitEuler": {
        "initialValues": iv,
        "numberTimeSteps": 5,
        "endTime": 0.1,
        "FiniteElementMethod": {
          "inputMeshIsGlobal": True,
          "nElements": [nx, ny],
          "physicalExtent": [2*nx, 2*ny],
          "relativeTolerance": 1e-15,
          "dirichletBoundaryConditions": {},
          "neumannBoundaryConditions": [],
        },
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out0", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings0(argc, argv, pythonConfig0);
  ProblemType problemSerial(settings0);
  problemSerial.run();


  std::vector<std::string> outputFilesToCheck = {"out0_0000004.py", "out0_0000004.0.py", "out0_0000004.1.py", "out0_0000004.2.py", "out0_0000004.3.py", "out0_0000004.4.py", "out0_0000004.5.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

// 2D regular fixed
TEST(DiffusionTest, SerialEqualsParallelLocal)
{
  typedef Control::MultipleInstances<
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<2>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<2>,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  > ProblemType;

  // run problem with global settings with only 1 rank
  std::string pythonConfig0 = R"(
# Diffusion 2D

import sys
rankNo = (int)(sys.argv[-2])
nRanks = (int)(sys.argv[-1])

print(rankNo,nRanks)

nx = 13   # number of elements in x direction
ny = 8   # number of elements in y direction

# 14*9=126 nodes

# initial values
iv = {}
iv[62] = 5.
iv[63] = 4.5
iv[64] = 4.
iv[79] = 4.5
iv[80] = 4.
iv[81] = 3.5

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [5],
      "ExplicitEuler": {
        "initialValues": iv,
        "numberTimeSteps": 5,
        "endTime": 0.1,
        "FiniteElementMethod": {
          "inputMeshIsGlobal": True,
          "nElements": [nx, ny],
          "physicalExtent": [2*nx, 2*ny],
          "relativeTolerance": 1e-15,
          "dirichletBoundaryConditions": {},
          "neumannBoundaryConditions": [],
          "prefactor": 1,
        },
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out1", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings0(argc, argv, pythonConfig0);
  int ownRankNo = settings0.ownRankNo();

  ProblemType problemSerial(settings0);
  problemSerial.run();

  // run problem with local settings
  std::string pythonConfig2 = R"(
# Diffusion 2D

import sys
rankNo = (int)(sys.argv[-2])
nRanks = (int)(sys.argv[-1])

print(rankNo,nRanks)

nx = 4   # number of elements in x direction
ny = 4   # number of elements in y direction

if rankNo == 0 or rankNo == 3:
  nx = 5

# 14*9=126 nodes

# initial values
iv = {}
if rankNo == 4:
  iv[1] = 5.  # 62
  iv[2] = 4.5 # 63
  iv[3] = 4.  # 64

elif rankNo == 5:
  iv[5] = 4.5  # 79
  iv[6] = 4.   # 80
  iv[7] = 3.5  # 81

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1,2,3,4,5],
      "ExplicitEuler": {
        "inputMeshIsGlobal": False,
        "initialValues": iv,
        "numberTimeSteps": 5,
        "endTime": 0.1,
        "FiniteElementMethod": {
          "inputMeshIsGlobal": False,
          "nElements": [nx, ny],
          "nRanks": [3,2],
          "physicalExtent": [2*nx, 2*ny],
          "relativeTolerance": 1e-15,
          "dirichletBoundaryConditions": {},
          "neumannBoundaryConditions": [],
        },
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out1", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  ProblemType problemLocalSettings(settings2);

  problemLocalSettings.run();

  std::vector<std::string> outputFilesToCheck = {"out1_0000004.py", "out1_0000004.0.py", "out1_0000004.1.py", "out1_0000004.2.py", "out1_0000004.3.py", "out1_0000004.4.py", "out1_0000004.5.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}
