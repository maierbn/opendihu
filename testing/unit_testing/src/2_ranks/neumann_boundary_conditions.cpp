#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"

namespace SpatialDiscretization
{

TEST(LaplaceTest, SerialEqualsParallelDeformable3DGlobal)
{
  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D, Neumann BC
import numpy as np
import sys

# 3D laplace problem

# run as (both with either local=True or local=False):
# mpirun -n 2 ./laplace ../settings_neumann.py
# ./laplace ../settings_neumann.py

local = False

nx = 2
ny = 2
nz = 5


# Neumann boundary conditions
bc = []

# global boundary conditions
if not local:
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      # z- plane (bottom)
      bc_value = np.sin(x*np.pi)
      bc_value = -1.0   # negative = inflow

      bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

      # z+ plane (top)
      element_no += int((nz-1)*nx*ny)
      bc_value = np.sin(y*np.pi) + 2.0
      bc_value = 1.0
      bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

rank_no = 0
n_ranks = 1

n_elements = [nx,ny,nz]

# local boundary conditions
if local:
  if n_ranks > 1:
    n_elements = [2,2,2]
    if rank_no == 0:
      n_elements = [2,2,3]

  # boundary conditions
  bc = []
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      if rank_no == 0:
        bc.append({"element": element_no, "constantValue": -1.0, "face": "2-"})

      if rank_no == n_ranks-1:
        bc.append({"element": -(nx*ny)+element_no, "constantValue": 1.0, "face": "2+"})

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod" : {
        "nElements": n_elements,
        "nRanks": [1,1,n_ranks],
        "inputMeshIsGlobal": not local,
        "physicalExtent": n_elements,
        "outputInterval": 1.0,
        "prefactor": 1,
        "dirichletBoundaryConditions": {0:0} if rank_no == 0 else {},
        "neumannBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "none",
        "maxIterations": 10000,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
        ]
      },
    }],
  }
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();

  LOG(INFO) << " =================== run parallel problem ================= ";

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 3D, Neumann BC
import numpy as np
import sys

# 3D laplace problem

# run as (both with either local=True or local=False):
# mpirun -n 2 ./laplace ../settings_neumann.py
# ./laplace ../settings_neumann.py

local = False

nx = 2
ny = 2
nz = 5


# Neumann boundary conditions
bc = []

# global boundary conditions
if not local:
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      # z- plane (bottom)
      bc_value = np.sin(x*np.pi)
      bc_value = -1.0   # negative = inflow

      bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

      # z+ plane (top)
      element_no += int((nz-1)*nx*ny)
      bc_value = np.sin(y*np.pi) + 2.0
      bc_value = 1.0
      bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

n_elements = [nx,ny,nz]

# local boundary conditions
if local:
  if n_ranks > 1:
    n_elements = [2,2,2]
    if rank_no == 0:
      n_elements = [2,2,3]

  # boundary conditions
  bc = []
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      if rank_no == 0:
        bc.append({"element": element_no, "constantValue": -1.0, "face": "2-"})

      if rank_no == n_ranks-1:
        bc.append({"element": -(nx*ny)+element_no, "constantValue": 1.0, "face": "2+"})

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod" : {
        "nElements": n_elements,
        "nRanks": [1,1,n_ranks],
        "inputMeshIsGlobal": not local,
        "physicalExtent": n_elements,
        "outputInterval": 1.0,
        "prefactor": 1,
        "dirichletBoundaryConditions": {0:0} if rank_no == 0 else {},
        "neumannBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "none",
        "maxIterations": 10000,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
        ]
      },
    }],
  }
}

)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out.py", "out.0.py", "out.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}
/*
TEST(LaplaceTest, SerialEqualsParallelDeformable3DLocal)
{
  // run serial problem
  std::string pythonConfig = R"(
# Laplace 3D, Neumann BC
import numpy as np
import sys

# 3D laplace problem

# run as (both with either local=True or local=False):
# mpirun -n 2 ./laplace ../settings_neumann.py
# ./laplace ../settings_neumann.py

local = True

nx = 2
ny = 2
nz = 5


# Neumann boundary conditions
bc = []

# global boundary conditions
if not local:
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      # z- plane (bottom)
      bc_value = np.sin(x*np.pi)
      bc_value = -1.0   # negative = inflow

      bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

      # z+ plane (top)
      element_no += int((nz-1)*nx*ny)
      bc_value = np.sin(y*np.pi) + 2.0
      bc_value = 1.0
      bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

rank_no = 0
n_ranks = 1

n_elements = [nx,ny,nz]

# local boundary conditions
if local:
  if n_ranks > 1:
    n_elements = [2,2,2]
    if rank_no == 0:
      n_elements = [2,2,3]

  # boundary conditions
  bc = []
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      if rank_no == 0:
        bc.append({"element": element_no, "constantValue": -1.0, "face": "2-"})

      if rank_no == n_ranks-1:
        bc.append({"element": -(nx*ny)+element_no, "constantValue": 1.0, "face": "2+"})

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0],
      "FiniteElementMethod" : {
        "nElements": n_elements,
        "nRanks": [1,1,n_ranks],
        "inputMeshIsGlobal": not local,
        "physicalExtent": n_elements,
        "outputInterval": 1.0,
        "prefactor": 1,
        "dirichletBoundaryConditions": {0:0} if rank_no == 0 else {},
        "neumannBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "none",
        "maxIterations": 10000,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
        ]
      },
    }],
  }
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  int ownRankNo = settings.ownRankNo();

  typedef Control::MultipleInstances<
    FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
    >
  > ProblemType;
  ProblemType problemSerial(settings);

  problemSerial.run();

  LOG(INFO) << " =================== run parallel problem ================= ";

  // run parallel problem
  std::string pythonConfig2 = R"(
# Laplace 3D, Neumann BC
import numpy as np
import sys

# 3D laplace problem

# run as (both with either local=True or local=False):
# mpirun -n 2 ./laplace ../settings_neumann.py
# ./laplace ../settings_neumann.py

local = True

nx = 2
ny = 2
nz = 5


# Neumann boundary conditions
bc = []

# global boundary conditions
if not local:
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      # z- plane (bottom)
      bc_value = np.sin(x*np.pi)
      bc_value = -1.0   # negative = inflow

      bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

      # z+ plane (top)
      element_no += int((nz-1)*nx*ny)
      bc_value = np.sin(y*np.pi) + 2.0
      bc_value = 1.0
      bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

n_elements = [nx,ny,nz]

# local boundary conditions
if local:
  if n_ranks > 1:
    n_elements = [2,2,2]
    if rank_no == 0:
      n_elements = [2,2,3]

  # boundary conditions
  bc = []
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)

      if rank_no == 0:
        bc.append({"element": element_no, "constantValue": -1.0, "face": "2-"})

      if rank_no == n_ranks-1:
        bc.append({"element": -(nx*ny)+element_no, "constantValue": 1.0, "face": "2+"})

config = {
  "MultipleInstances": {
    "nInstances": 1,
    "instances": [{
      "ranks": [0,1],
      "FiniteElementMethod" : {
        "nElements": n_elements,
        "nRanks": [1,1,n_ranks],
        "inputMeshIsGlobal": not local,
        "physicalExtent": n_elements,
        "outputInterval": 1.0,
        "prefactor": 1,
        "dirichletBoundaryConditions": {0:0} if rank_no == 0 else {},
        "neumannBoundaryConditions": bc,
        "relativeTolerance": 1e-15,
        "solverType": "gmres",
        "preconditionerType": "none",
        "maxIterations": 10000,
        "OutputWriter" : [
          {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
        ]
      },
    }],
  }
}

)";


  DihuContext settings2(argc, argv, pythonConfig2);

  ProblemType problemParallel(settings2);

  problemParallel.run();

  std::vector<std::string> outputFilesToCheck = {"out.py", "out.0.py", "out.1.py"};
  if (ownRankNo == 0)
  {
    assertParallelEqualsSerialOutputFiles(outputFilesToCheck);
  }

  nFails += ::testing::Test::HasFailure();
}
*/
}  // namespace

