#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"

// the following tests only fail on travis ci but succeed anywhere else
#ifndef ON_TRAVIS_CI

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
