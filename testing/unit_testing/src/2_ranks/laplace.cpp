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
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
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

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000004, 0.7999999999999998, 0.6000000000000001]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -1.25, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [3], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.4000000000000001, 0.1999999999999999, 0.0]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  assertFileMatchesContent("out.0.py", referenceOutput0);
  assertFileMatchesContent("out.1.py", referenceOutput1);

  nFails += ::testing::Test::HasFailure();
}

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
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
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

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Lagrange\", \"basisOrder\": 2, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.4, 0.8, 1.2000000000000002, 1.6, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000004, 0.8999999999999994, 0.799999999999999, 0.699999999999998, 0.5999999999999972, 0.49999999999999767]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -3.333333333333333, 0.41666666666666674, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [6], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 2, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 2.8000000000000003, 3.2, 3.6, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.39999999999999736, 0.2999999999999976, 0.19999999999999823, 0.09999999999999895, 0.0]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  assertFileMatchesContent("out.0.py", referenceOutput0);
  assertFileMatchesContent("out.1.py", referenceOutput1);

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
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False, "onlyNodalValues": True}
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

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000002, 0.8000000000000003, 0.5999999999999999]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -1.4999999999999993, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [3], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.3999999999999998, 0.1999999999999998, 0.0]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  assertFileMatchesContent("out.0.py", referenceOutput0);
  assertFileMatchesContent("out.1.py", referenceOutput1);

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
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "filename": "out", "outputInterval": 1, "binary": False},
      {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
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

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [3], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [false], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.02933938631260563, 0.0046848256044837895, -0.0014240802367548704]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [-0.031249999999999965, 0.031249999999999965, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [2], \"beginNodeGlobalNatural\": [3], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Hermite\", \"basisOrder\": 3, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.00039567686177667104, 6.594614362942627e-05, 0.0]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  assertFileMatchesContent("out.0.py", referenceOutput0);
  assertFileMatchesContent("out.1.py", referenceOutput1);

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
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
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

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [3, 2], \"nElementsLocal\": [2, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [false, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 2.0, 0.0, 2.0, 0.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 2.0, 2.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0000000000000002, 4.242857142857144, 5.971428571428574, 0.0, 10.000000000000005]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0, -3.666666666666667, -11.000000000000002, 0.0, 10.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput1 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [3, 2], \"nElementsLocal\": [1, 2], \"beginNodeGlobalNatural\": [2, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 2, \"ownRankNo\": 1, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [4.0, 6.0, 4.0, 6.0, 4.0, 6.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 2.0, 2.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [2.0000000000000004, 3.0000000000000018, 10.528571428571434, 12.257142857142862, 20.00000000000001, 30.000000000000007]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [2.0, 3.0, -22.0, -12.833333333333332, 20.0, 30.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  assertFileMatchesContent("out2d_p2.0.py", referenceOutput0);
  assertFileMatchesContent("out2d_p2.1.py", referenceOutput1);

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
        "DirichletBoundaryCondition": bc,
        "relativeTolerance": 1e-15,
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

  std::string referenceOutput0 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [3, 2], \"nElementsLocal\": [3, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0, 2.0, 3.0, 4.242857142857142, 5.971428571428571, 10.52857142857143, 12.257142857142851, 0.0, 9.999999999999996, 19.999999999999993, 30.000000000000004]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [0.0, 1.0, 2.0, 3.0, -3.666666666666667, -11.000000000000002, -22.0, -12.833333333333332, 0.0, 10.0, 20.0, 30.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";

  assertFileMatchesContent("out2d_p1.py", referenceOutput0);

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, Structured2DLinearParallelWithMultipleInstances)
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
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [6.0, 4.0],
        "DirichletBoundaryCondition": bc,
        "relativeTolerance": 1e-15,
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

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelStructured2DLinear)
{
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

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out.py", "out.0.py", "out.1.py"};
  assertParallelEqualsSerialOutputFiles(outputFilesToCheck);

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelStructured2DQuadratic)
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

TEST(LaplaceTest, SerialEqualsParallelStructured2DHermite)
{
  // run serial problem
  std::string pythonConfig = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(2*nx+1)):
  x = i/(2*nx+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*nx+1)*ny + i
  bc[i2] = np.sin(x*np.pi)

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
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
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<2>,
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

# boundary conditions
bc = {}
for i in range(int(2*nx+1)):
  x = i/(2*nx+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*nx+1)*ny + i
  bc[i2] = np.sin(x*np.pi)

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [2*nx, 2*ny],
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

// structured deformable

TEST(LaplaceTest, SerialEqualsParallelRegular2DLinear)
{
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

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
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

  Control::MultipleInstances<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Static::Laplace
    >
  > problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out.py", "out.0.py", "out.1.py"};
  assertParallelEqualsSerialOutputFiles(outputFilesToCheck);

  nFails += ::testing::Test::HasFailure();
}

TEST(LaplaceTest, SerialEqualsParallelRegular2DQuadratic)
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

TEST(LaplaceTest, SerialEqualsParallelRegular2DHermite)
{
  // run serial problem
  std::string pythonConfig = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(2*nx+1)):
  x = i/(2*nx+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*nx+1)*ny + i
  bc[i2] = np.sin(x*np.pi)

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

  // run parallel problem
  std::string pythonConfig2 = R"(
import numpy as np

nx = 10   # number of elements in x direction
ny = 12   # number of elements in y direction

# boundary conditions
bc = {}
for i in range(int(2*nx+1)):
  x = i/(2*nx+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*nx+1)*ny + i
  bc[i2] = np.sin(x*np.pi)

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
