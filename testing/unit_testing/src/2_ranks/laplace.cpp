#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"

TEST(LaplaceTest, Structured1DLinear)
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
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out4", "outputInterval": 1, "binary": False, "fixedFormat": True, "combineFiles": False, "onlyNodalValues": True},
      {"format": "PythonFile", "filename": "out4", "outputInterval": 1, "binary": False, "onlyNodalValues": True}
    ]
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

  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000004, 0.7999999999999998, 0.6000000000000001]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -1.25, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [3], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.4000000000000001, 0.1999999999999999, 0.0]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  if (settings.ownRankNo() == 0)
  {
    assertFileMatchesContent("out4.0.py", referenceOutput0);
    assertFileMatchesContent("out4.1.py", referenceOutput1);
  }

  nFails += ::testing::Test::HasFailure();
}
/*
TEST(LaplaceTest, Structured1DQuadratic)
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
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out5", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out5", "outputInterval": 1, "binary": False}
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > problem(settings);

  problem.run();

  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Lagrange\", \"basisOrder\": 2, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.4, 0.8, 1.2000000000000002, 1.6, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000004, 0.8999999999999994, 0.799999999999999, 0.699999999999998, 0.5999999999999972, 0.49999999999999767]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -3.333333333333333, 0.41666666666666674, 0.0, 0.0, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [6], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 2, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 2.8000000000000003, 3.2, 3.6, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.39999999999999736, 0.2999999999999976, 0.19999999999999823, 0.09999999999999895, 0.0]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  if (settings.ownRankNo() == 0)
  {
    assertFileMatchesContent("out5.0.py", referenceOutput0);
    assertFileMatchesContent("out5.1.py", referenceOutput1);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, Structured1DHermite1)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5

# boundary conditions
bc = {}
bc[0] = 1.0
bc[-2] = 0.0

config = {
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out6", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out6", "outputInterval": 1, "binary": False, "onlyNodalValues": True}
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > problem(settings);

  problem.run();

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000002, 0.8000000000000003, 0.5999999999999999]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -1.4999999999999993, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [3], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.3999999999999998, 0.1999999999999998, 0.0]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  if (settings.ownRankNo() == 0)
  {
    assertFileMatchesContent("out6.0.py", referenceOutput0);
    assertFileMatchesContent("out6.1.py", referenceOutput1);
  }

  nFails += ::testing::Test::HasFailure();

}

TEST(LaplaceTest, Structured1DHermite2)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5

# boundary conditions
bc = {}
bc[1] = -1/4.
bc[-2] = 0.0

config = {
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out7", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out7", "outputInterval": 1, "binary": False, "onlyNodalValues": False}
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > problem(settings);

  problem.run();

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": false, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.0, 0.8, 0.0, 1.6, 0.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [-1.731137919321699e-17, -0.24999999999999994, -2.2431623986420583e-18, -0.05628740946730341, 1.3321963259718046e-18, -0.01267457751408284]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [1.3877787807814457e-17, -0.25, -1.3877787807814457e-17, 0.026041666666666668, 0.0, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [3], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": false, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 0.0, 3.2, 0.0, 4.0, 0.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.2293431050318027e-18, -0.002860618931749887, 1.0301064917240598e-18, -0.000674977500749973, 0.0, -0.0002892760717499882]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  if (settings.ownRankNo() == 0)
  {
    assertFileMatchesContent("out7.0.py", referenceOutput0);
    assertFileMatchesContent("out7.1.py", referenceOutput1);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, Structured2DLinear)
{
  std::string pythonConfig = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "nElements": [nx, ny],
    "physicalExtent": [6.0, 4.0],
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out2d_p2", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out2d_p2", "outputInterval": 1, "binary": False}
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > problem(settings);

  problem.run();

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [3, 2], \"nElementsLocal\": [2, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [false, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 2.0, 0.0, 2.0, 0.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 2.0, 2.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0000000000000002, 4.242857142857144, 5.971428571428574, 0.0, 10.000000000000005]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0, -3.666666666666667, -11.000000000000002, 0.0, 10.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [3, 2], \"nElementsLocal\": [1, 2], \"beginNodeGlobalNatural\": [2, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [4.0, 6.0, 4.0, 6.0, 4.0, 6.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 2.0, 2.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [2.0000000000000004, 3.0000000000000018, 10.528571428571434, 12.257142857142862, 20.00000000000001, 30.000000000000007]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [2.0, 3.0, -22.0, -12.833333333333332, 20.0, 30.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  if (settings.ownRankNo() == 0)
  {
    assertFileMatchesContent("out2d_p2.0.py", referenceOutput0);
    assertFileMatchesContent("out2d_p2.1.py", referenceOutput1);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, Structured2DLinearSerial)
{
  std::string pythonConfig = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "Paraview", "filename": "out2d_p1", "outputInterval": 1, "binary": False},
          {"format": "PythonFile", "filename": "out2d_p1", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [3, 2], \"nElementsLocal\": [3, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0, 2.0, 3.0, 4.242857142857142, 5.971428571428571, 10.52857142857143, 12.257142857142851, 0.0, 9.999999999999996, 19.999999999999993, 30.000000000000004]}]}, {\"name\": \"rightHandSide\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0, 2.0, 3.0, -3.666666666666667, -11.000000000000002, -22.0, -12.833333333333332, 0.0, 10.0, 20.0, 30.0]}]}, {\"name\": \"-rhsNeumannBC\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  if (settings.ownRankNo() == 0)
  {
    assertFileMatchesContent("out2d_p1.py", referenceOutput0);
  }

  nFails += ::testing::Test::HasFailure();
}

// the following tests only fail on travis ci but succeed anywhere else
#ifndef ON_TRAVIS_CI
TEST(LaplaceTest, Structured2DLinearParallelWithMultipleInstances)
{
  std::cout << "wait 1 s" << std::endl;
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  std::string pythonConfig = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "Paraview", "filename": "out2d_p1", "outputInterval": 1, "binary": False},
          {"format": "PythonFile", "filename": "out2d_p1", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problem(settings);
  LOG(INFO) << "problem created";

  //problem.run();

  nFails += ::testing::Test::HasFailure();
}


// 2D structured deformable
TEST(LaplaceTest, SerialEqualsParallelRegular2DLinear)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out8", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemSerial(settings);

  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "Paraview", "filename": "out8", "outputInterval": 1, "binary": False},
          {"format": "PythonFile", "filename": "out8", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out8.py", "out8.0.py", "out8.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelRegular2DQuadratic)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
import numpy as np

nx = 5   # number of elements in x direction
ny = 6   # number of elements in y direction

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
        "physicalExtent": [2*nx, 2*ny],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out9", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
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

nx = 5   # number of elements in x direction
ny = 6   # number of elements in y direction

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
        "physicalExtent": [2*nx, 2*ny],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out9", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out9.py", "out9.0.py", "out9.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelRegular2DHermite)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny

n_nodes_x = nx+1
n_nodes_y = ny+1

# f(x,y) = x^2 - y^2
# grad f = [2*x; -2*y]
# f_xx = 2
# f_yy = -2
# Δf = 2-2 = 0
#
# boundary value
def f(x,y):
  return x**2 - y**2

def f_x(x,y):
  return 2*x

def f_y(x,y):
  return -2*y

# f_xy is zero

# boundary conditions
bc = {}
# y-, y+
for i in range(n_nodes_x):
  x = float(i)/(n_nodes_x-1) * physical_extend_x

  # y-
  y = 0
  bc[4*(i)+0] = f(x,y)
  bc[4*(i)+1] = f_x(x,y)
  bc[4*(i)+2] = f_y(x,y)
  bc[4*(i)+3] = 0   # f_xy

  # y+
  y = physical_extend_y
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+0] = f(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+1] = f_x(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+2] = f_y(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+3] = 0   # f_xy

# x-, x+
for j in range(n_nodes_y):
  y = float(j)/(n_nodes_y-1) * physical_extend_y

  # x-
  x = 0
  bc[4*(j*n_nodes_x)+0] = f(x,y)
  bc[4*(j*n_nodes_x)+1] = f_x(x,y)
  bc[4*(j*n_nodes_x)+2] = f_y(x,y)
  bc[4*(j*n_nodes_x)+3] = 0   # f_xy

  # x+
  x = physical_extend_x
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+0] = f(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+1] = f_x(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+2] = f_y(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+3] = 0   # f_xy

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out10", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::Hermite,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny

n_nodes_x = nx+1
n_nodes_y = ny+1

# f(x,y) = x^2 - y^2
# grad f = [2*x; -2*y]
# f_xx = 2
# f_yy = -2
# Δf = 2-2 = 0
#
# boundary value
def f(x,y):
  return x**2 - y**2

def f_x(x,y):
  return 2*x

def f_y(x,y):
  return -2*y

# f_xy is zero

# boundary conditions
bc = {}
# y-, y+
for i in range(n_nodes_x):
  x = float(i)/(n_nodes_x-1) * physical_extend_x

  # y-
  y = 0
  bc[4*(i)+0] = f(x,y)
  bc[4*(i)+1] = f_x(x,y)
  bc[4*(i)+2] = f_y(x,y)
  bc[4*(i)+3] = 0   # f_xy

  # y+
  y = physical_extend_y
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+0] = f(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+1] = f_x(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+2] = f_y(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+3] = 0   # f_xy

# x-, x+
for j in range(n_nodes_y):
  y = float(j)/(n_nodes_y-1) * physical_extend_y

  # x-
  x = 0
  bc[4*(j*n_nodes_x)+0] = f(x,y)
  bc[4*(j*n_nodes_x)+1] = f_x(x,y)
  bc[4*(j*n_nodes_x)+2] = f_y(x,y)
  bc[4*(j*n_nodes_x)+3] = 0   # f_xy

  # x+
  x = physical_extend_x
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+0] = f(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+1] = f_x(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+2] = f_y(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+3] = 0   # f_xy

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out10", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out10.py", "out10.0.py", "out10.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }
  else
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));  // pause execution, such that output files can be closed
  }

  nFails += ::testing::Test::HasFailure();
}

// 2D structured regular fixed
TEST(LaplaceTest, SerialEqualsParallelDeformable2DLinear)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out11", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemSerial(settings);

  problemSerial.run();

  std::this_thread::sleep_for(std::chrono::milliseconds(100));  // pause execution, such that output files can be closed

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 2D, 3 x 2 (=6) elements, 4 x 3 (=12) nodes

import numpy as np

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = i
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 10*i

print("{} x {} nodes, Dirichlet boundary conditions: {}".format(nx+1,ny+1,bc))

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out11", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out11.py", "out11.0.py", "out11.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }
  else
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));  // pause execution, such that output files can be closed
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelDeformable2DQuadratic)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

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
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out12", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";
  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

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
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out12", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out12.py", "out12.0.py", "out12.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }
  else
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));  // pause execution, such that output files can be closed
  }

  nFails += ::testing::Test::HasFailure();
}


TEST(LaplaceTest, SerialEqualsParallelDeformable2DHermite)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

physical_extend_x = 6.0
physical_extend_y = 4.0

n_nodes_x = nx+1
n_nodes_y = ny+1

# f(x,y) = x^2 - y^2
# grad f = [2*x; -2*y]
# f_xx = 2
# f_yy = -2
# Δf = 2-2 = 0
#
# boundary value
def f(x,y):
  return x**2 - y**2

def f_x(x,y):
  return 2*x

def f_y(x,y):
  return -2*y

# f_xy is zero

# boundary conditions
bc = {}
# y-, y+
for i in range(n_nodes_x):
  x = float(i)/(n_nodes_x-1) * physical_extend_x

  # y-
  y = 0
  bc[4*(i)+0] = f(x,y)
  bc[4*(i)+1] = f_x(x,y)
  bc[4*(i)+2] = f_y(x,y)
  bc[4*(i)+3] = 0   # f_xy

  # y+
  y = physical_extend_y
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+0] = f(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+1] = f_x(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+2] = f_y(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+3] = 0   # f_xy

# x-, x+
for j in range(n_nodes_y):
  y = float(j)/(n_nodes_y-1) * physical_extend_y

  # x-
  x = 0
  bc[4*(j*n_nodes_x)+0] = f(x,y)
  bc[4*(j*n_nodes_x)+1] = f_x(x,y)
  bc[4*(j*n_nodes_x)+2] = f_y(x,y)
  bc[4*(j*n_nodes_x)+3] = 0   # f_xy

  # x+
  x = physical_extend_x
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+0] = f(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+1] = f_x(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+2] = f_y(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+3] = 0   # f_xy

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [physical_extend_x, physical_extend_y],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out13", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::Hermite,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();


  std::this_thread::sleep_for(std::chrono::milliseconds(100));  // pause execution, such that output files can be closed

  // run parallel problem
  std::string pythonConfig2 = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

physical_extend_x = 6.0
physical_extend_y = 4.0

n_nodes_x = nx+1
n_nodes_y = ny+1

# f(x,y) = x^2 - y^2
# grad f = [2*x; -2*y]
# f_xx = 2
# f_yy = -2
# Δf = 2-2 = 0
#
# boundary value
def f(x,y):
  return x**2 - y**2

def f_x(x,y):
  return 2*x

def f_y(x,y):
  return -2*y

# f_xy is zero

# boundary conditions
bc = {}
# y-, y+
for i in range(n_nodes_x):
  x = float(i)/(n_nodes_x-1) * physical_extend_x

  # y-
  y = 0
  bc[4*(i)+0] = f(x,y)
  bc[4*(i)+1] = f_x(x,y)
  bc[4*(i)+2] = f_y(x,y)
  bc[4*(i)+3] = 0   # f_xy

  # y+
  y = physical_extend_y
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+0] = f(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+1] = f_x(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+2] = f_y(x,y)
  bc[4*((n_nodes_y-1)*n_nodes_x + i)+3] = 0   # f_xy

# x-, x+
for j in range(n_nodes_y):
  y = float(j)/(n_nodes_y-1) * physical_extend_y

  # x-
  x = 0
  bc[4*(j*n_nodes_x)+0] = f(x,y)
  bc[4*(j*n_nodes_x)+1] = f_x(x,y)
  bc[4*(j*n_nodes_x)+2] = f_y(x,y)
  bc[4*(j*n_nodes_x)+3] = 0   # f_xy

  # x+
  x = physical_extend_x
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+0] = f(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+1] = f_x(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+2] = f_y(x,y)
  bc[4*(j*n_nodes_x + (n_nodes_x-1))+3] = 0   # f_xy

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [physical_extend_x, physical_extend_y],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out13", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out13.py", "out13.0.py", "out13.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

// 3D structured deformable
TEST(LaplaceTest, SerialEqualsParallelDeformable3DLinear)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[j*(nx+1)+i] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[i2] = 10.*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out14", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemSerial(settings);

  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[j*(nx+1)+i] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[i2] = 10.*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out14", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out14.py", "out14.0.py", "out14.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }
  else
  {

  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelDeformable3DQuadratic)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny
physical_extend_z = 2*nz

n_dofs_x = 2*nx+1
n_dofs_y = 2*ny+1
n_dofs_z = 2*nz+1

# f(x,y,z) = x^2 - 1/2*y^2 - 1/2*z^2
# grad f = [2*x; -y; -z]
# Δf = 2-1-1 = 0
#
# boundary value
def solution(x,y,z):
  return x**2 - 0.5*y**2 - 0.5*z**2


# boundary conditions
bc = {}
# z-, z+
for j in range(n_dofs_y):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # z-
    z = 0
    bc[j*n_dofs_x + i] = solution(x,y,z)

    # z+
    z = physical_extend_z
    bc[(n_dofs_z-1)*n_dofs_y*n_dofs_x + j*n_dofs_x + i] = solution(x,y,z)

# y-, y+
for k in range(n_dofs_z):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    z = float(k)/(n_dofs_z-1) * physical_extend_z

    # y-
    y = 0
    bc[k*n_dofs_y*n_dofs_x + i] = solution(x,y,z)

    # y+
    y = physical_extend_y
    bc[k*n_dofs_y*n_dofs_x + (n_dofs_y-1)*n_dofs_x + i] = solution(x,y,z)

# x-, x+
for k in range(n_dofs_z):
  for j in range(n_dofs_y):
    z = float(k)/(n_dofs_z-1) * physical_extend_z
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # x-
    x = 0
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x] = solution(x,y,z)

    # x+
    x = physical_extend_x
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x + (n_dofs_x-1)] = solution(x,y,z)

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [physical_extend_x, physical_extend_y, physical_extend_z],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out15", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<3>,
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
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny
physical_extend_z = 2*nz

n_dofs_x = 2*nx+1
n_dofs_y = 2*ny+1
n_dofs_z = 2*nz+1

# f(x,y,z) = x^2 - 1/2*y^2 - 1/2*z^2
# grad f = [2*x; -y; -z]
# Δf = 2-1-1 = 0
#
# boundary value
def solution(x,y,z):
  return x**2 - 0.5*y**2 - 0.5*z**2


# boundary conditions
bc = {}
# z-, z+
for j in range(n_dofs_y):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # z-
    z = 0
    bc[j*n_dofs_x + i] = solution(x,y,z)

    # z+
    z = physical_extend_z
    bc[(n_dofs_z-1)*n_dofs_y*n_dofs_x + j*n_dofs_x + i] = solution(x,y,z)

# y-, y+
for k in range(n_dofs_z):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    z = float(k)/(n_dofs_z-1) * physical_extend_z

    # y-
    y = 0
    bc[k*n_dofs_y*n_dofs_x + i] = solution(x,y,z)

    # y+
    y = physical_extend_y
    bc[k*n_dofs_y*n_dofs_x + (n_dofs_y-1)*n_dofs_x + i] = solution(x,y,z)

# x-, x+
for k in range(n_dofs_z):
  for j in range(n_dofs_y):
    z = float(k)/(n_dofs_z-1) * physical_extend_z
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # x-
    x = 0
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x] = solution(x,y,z)

    # x+
    x = physical_extend_x
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x + (n_dofs_x-1)] = solution(x,y,z)

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out15", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out15.py", "out15.0.py", "out15.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

// Test does not converge and gives slightly different results
TEST(LaplaceTest, SerialEqualsParallelDeformable3DHermite)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[8*(j*(nx+1)+i)] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[8*i2] = 10.*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "maxIterations": 1e6,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out16", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

 int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::Hermite,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[8*(j*(nx+1)+i)] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[8*i2] = 10.*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "maxIterations": 1e5,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out16", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out16.py", "out16.0.py", "out16.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

// 3D structured regular fixed
TEST(LaplaceTest, SerialEqualsParallelRegular3DLinear)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[j*(nx+1)+i] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[i2] = 10.*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out17", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemSerial(settings);

  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[j*(nx+1)+i] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[i2] = 10.*i

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out17", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out17.py", "out17.0.py", "out17.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelRegular3DQuadratic)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny
physical_extend_z = 2*nz

n_dofs_x = 2*nx+1
n_dofs_y = 2*ny+1
n_dofs_z = 2*nz+1

# f(x,y,z) = x^2 - 1/2*y^2 - 1/2*z^2
# grad f = [2*x; -y; -z]
# Δf = 2-1-1 = 0
#
# boundary value
def solution(x,y,z):
  return x**2 - 0.5*y**2 - 0.5*z**2


# boundary conditions
bc = {}
# z-, z+
for j in range(n_dofs_y):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # z-
    z = 0
    bc[j*n_dofs_x + i] = solution(x,y,z)

    # z+
    z = physical_extend_z
    bc[(n_dofs_z-1)*n_dofs_y*n_dofs_x + j*n_dofs_x + i] = solution(x,y,z)

# y-, y+
for k in range(n_dofs_z):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    z = float(k)/(n_dofs_z-1) * physical_extend_z

    # y-
    y = 0
    bc[k*n_dofs_y*n_dofs_x + i] = solution(x,y,z)

    # y+
    y = physical_extend_y
    bc[k*n_dofs_y*n_dofs_x + (n_dofs_y-1)*n_dofs_x + i] = solution(x,y,z)

# x-, x+
for k in range(n_dofs_z):
  for j in range(n_dofs_y):
    z = float(k)/(n_dofs_z-1) * physical_extend_z
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # x-
    x = 0
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x] = solution(x,y,z)

    # x+
    x = physical_extend_x
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x + (n_dofs_x-1)] = solution(x,y,z)

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [physical_extend_x, physical_extend_y, physical_extend_z],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out18", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<3>,
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
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny
physical_extend_z = 2*nz

n_dofs_x = 2*nx+1
n_dofs_y = 2*ny+1
n_dofs_z = 2*nz+1

# f(x,y,z) = x^2 - 1/2*y^2 - 1/2*z^2
# grad f = [2*x; -y; -z]
# Δf = 2-1-1 = 0
#
# boundary value
def solution(x,y,z):
  return x**2 - 0.5*y**2 - 0.5*z**2


# boundary conditions
bc = {}
# z-, z+
for j in range(n_dofs_y):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # z-
    z = 0
    bc[j*n_dofs_x + i] = solution(x,y,z)

    # z+
    z = physical_extend_z
    bc[(n_dofs_z-1)*n_dofs_y*n_dofs_x + j*n_dofs_x + i] = solution(x,y,z)

# y-, y+
for k in range(n_dofs_z):
  for i in range(n_dofs_x):
    x = float(i)/(n_dofs_x-1) * physical_extend_x
    z = float(k)/(n_dofs_z-1) * physical_extend_z

    # y-
    y = 0
    bc[k*n_dofs_y*n_dofs_x + i] = solution(x,y,z)

    # y+
    y = physical_extend_y
    bc[k*n_dofs_y*n_dofs_x + (n_dofs_y-1)*n_dofs_x + i] = solution(x,y,z)

# x-, x+
for k in range(n_dofs_z):
  for j in range(n_dofs_y):
    z = float(k)/(n_dofs_z-1) * physical_extend_z
    y = float(j)/(n_dofs_y-1) * physical_extend_y

    # x-
    x = 0
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x] = solution(x,y,z)

    # x+
    x = physical_extend_x
    bc[k*n_dofs_y*n_dofs_x + j*n_dofs_x + (n_dofs_x-1)] = solution(x,y,z)
config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out18", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out18.py", "out18.0.py", "out18.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}


// Test does sometimes not converge and gives slightly different solutions
TEST(LaplaceTest, SerialEqualsParallelRegular3DHermite)
{
  LOG(INFO) << "wait 1 s";
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));  // pause execution, such that output files can be closed

  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny
physical_extend_z = 2*nz

n_nodes_x = nx+1
n_nodes_y = ny+1
n_nodes_z = nz+1

# f(x,y,z) = x^2 - 1/2*y^2 - 1/2*z^2
# grad f = [2*x; -y; -z]
# f_xx = 2
# f_xy
# Δf = 2-1-1 = 0
#
# boundary value
def f(x,y,z):
  return x**2 - 0.5*y**2 - 0.5*z**2

def f_x(x,y,z):
  return 2*x

def f_y(x,y,z):
  return -y

def f_z(x,y,z):
  return -z

# f_xy, f_xz, f_yz, f_xyz are zero

# boundary conditions
bc = {}
# z-, z+
for j in range(n_nodes_y):
  for i in range(n_nodes_x):
    x = float(i)/(n_nodes_x-1) * physical_extend_x
    y = float(j)/(n_nodes_y-1) * physical_extend_y

    # z-
    z = 0
    bc[8*(j*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*(j*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*(j*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*(j*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*(j*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*(j*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*(j*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*(j*n_nodes_x + i)+7] = 0   # f_zyx

    # z+
    z = physical_extend_z
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+7] = 0   # f_zyx

# y-, y+
for k in range(n_nodes_z):
  for i in range(n_nodes_x):
    x = float(i)/(n_nodes_x-1) * physical_extend_x
    z = float(k)/(n_nodes_z-1) * physical_extend_z

    # y-
    y = 0
    bc[8*(k*n_nodes_y*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + i)+7] = 0   # f_zyx

    # y+
    y = physical_extend_y
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+7] = 0   # f_zyx

# x-, x+
for k in range(n_nodes_z):
  for j in range(n_nodes_y):
    z = float(k)/(n_nodes_z-1) * physical_extend_z
    y = float(j)/(n_nodes_y-1) * physical_extend_y

    # x-
    x = 0
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+7] = 0   # f_zyx

    # x+
    x = physical_extend_x
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+7] = 0   # f_zyx

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "maxIterations": 1e5,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out19", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

 int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<3>,
      BasisFunction::Hermite,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 3D

nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

physical_extend_x = 2*nx
physical_extend_y = 2*ny
physical_extend_z = 2*nz

n_nodes_x = nx+1
n_nodes_y = ny+1
n_nodes_z = nz+1

# f(x,y,z) = x^2 - 1/2*y^2 - 1/2*z^2
# grad f = [2*x; -y; -z]
# f_xx = 2
# f_xy
# Δf = 2-1-1 = 0
#
# boundary value
def f(x,y,z):
  return x**2 - 0.5*y**2 - 0.5*z**2

def f_x(x,y,z):
  return 2*x

def f_y(x,y,z):
  return -y

def f_z(x,y,z):
  return -z

# f_xy, f_xz, f_yz, f_xyz are zero

# boundary conditions
bc = {}
# z-, z+
for j in range(n_nodes_y):
  for i in range(n_nodes_x):
    x = float(i)/(n_nodes_x-1) * physical_extend_x
    y = float(j)/(n_nodes_y-1) * physical_extend_y

    # z-
    z = 0
    bc[8*(j*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*(j*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*(j*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*(j*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*(j*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*(j*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*(j*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*(j*n_nodes_x + i)+7] = 0   # f_zyx

    # z+
    z = physical_extend_z
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*((n_nodes_z-1)*n_nodes_y*n_nodes_x + j*n_nodes_x + i)+7] = 0   # f_zyx

# y-, y+
for k in range(n_nodes_z):
  for i in range(n_nodes_x):
    x = float(i)/(n_nodes_x-1) * physical_extend_x
    z = float(k)/(n_nodes_z-1) * physical_extend_z

    # y-
    y = 0
    bc[8*(k*n_nodes_y*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + i)+7] = 0   # f_zyx

    # y+
    y = physical_extend_y
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + (n_nodes_y-1)*n_nodes_x + i)+7] = 0   # f_zyx

# x-, x+
for k in range(n_nodes_z):
  for j in range(n_nodes_y):
    z = float(k)/(n_nodes_z-1) * physical_extend_z
    y = float(j)/(n_nodes_y-1) * physical_extend_y

    # x-
    x = 0
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x)+7] = 0   # f_zyx

    # x+
    x = physical_extend_x
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+0] = f(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+1] = f_x(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+2] = f_y(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+3] = 0   # f_xy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+4] = f_z(x,y,z)
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+5] = 0   # f_zx
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+6] = 0   # f_zy
    bc[8*(k*n_nodes_y*n_nodes_x + j*n_nodes_x + (n_nodes_x-1))+7] = 0   # f_zyx

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [2*nx, 2*ny, 2*nz],
        "dirichletBoundaryConditions": bc,
        "maxIterations": 1e5,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "sor",
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out19", "outputInterval": 1, "binary": False}
        ]
      }
    }]
  }
}
)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out19.py", "out19.0.py", "out19.1.py"};
  if (ownRankNo == 0)
 {
   assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
 }

  nFails += ::testing::Test::HasFailure();
}
#endif
// */
