#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

namespace SpatialDiscretization
{

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect1DLinear)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5

# Neumann boundary conditions
bc = [
  {"element": 0, "constantValue": -1.0, "face": "0-"},
  {"element": n-1, "constantValue": 1.0, "face": "0+"}
]
config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 1.0,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "neumannBoundaryConditions": bc,
    "dirichletBoundaryConditions": {1:0},
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
    ]
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // ---- linear -----
  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized1(settings);
  equationDiscretized1.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized2(settings);
  equationDiscretized2.run();

  // unstructured
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized3(settings);
  equationDiscretized3.run();

  std::vector<double> referenceSolutionLinear = {
    -0.2, 0, 0.2, 0.4, 0.6, 0.8
  };
  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionLinear);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionLinear);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionLinear);
}


TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionBackwardsIsCorrect1DLinear)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5

# Neumann boundary conditions
bc = [
  {"element": 0, "constantValue": -1.0, "face": "0-"},
  {"element": -1, "constantValue": 1.0, "face": "0+"}
]
config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 1.0,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "neumannBoundaryConditions": bc,
    "dirichletBoundaryConditions": {1:0},
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
    ]
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // ---- linear -----
  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized1(settings);
  equationDiscretized1.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized2(settings);
  equationDiscretized2.run();

  std::vector<double> referenceSolutionLinear = {
    -0.2, 0, 0.2, 0.4, 0.6, 0.8
  };
  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionLinear);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionLinear);
}


TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect1DHermite)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5

# Neumann boundary conditions
bc = [
  {"element": 0, "constantValue": -1.0, "face": "0-"},
  {"element": n-1, "constantValue": 1.0, "face": "0+"}
]
config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 1.0,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "neumannBoundaryConditions": bc,
    "dirichletBoundaryConditions": {0:0, 1:1./5},
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
    ]
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // ---- Hermite -----
  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::Hermite,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized4(settings);
  equationDiscretized4.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::Hermite,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized5(settings);
  equationDiscretized5.run();

  // unstructured
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<1>,
    BasisFunction::Hermite,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized6(settings);
  equationDiscretized6.run();

  std::vector<double> referenceSolutionHermite = {
     0,  0.2,  0.2,  0.2,  0.4,  0.2,  0.6,  0.2,  0.8,  0.2,  1,  0.2
  };
  StiffnessMatrixTester::compareSolution(equationDiscretized4, referenceSolutionHermite);
  StiffnessMatrixTester::compareSolution(equationDiscretized5, referenceSolutionHermite);
  StiffnessMatrixTester::compareSolution(equationDiscretized6, referenceSolutionHermite);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect1DQuadratic)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5

# Neumann boundary conditions
bc = [
  {"element": 0, "constantValue": -1.0, "face": "0-"},
  {"element": n-1, "constantValue": 1.0, "face": "0+"}
]
config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 1.0,
    "nodePositions": [[float(i)/(2*n)] for i in range(2*n+1)],
    "elements": [[2*i, 2*i+1, 2*i+2] for i in range(n)],
    "neumannBoundaryConditions": bc,
    "dirichletBoundaryConditions": {1:0},
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
    ]
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // ---- quadratic -----
  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized4(settings);
  equationDiscretized4.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized5(settings);
  equationDiscretized5.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized6(settings);
  equationDiscretized6.run();

  std::vector<double> referenceSolutionQuadratic = {
    -0.1, 0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9
  };
  StiffnessMatrixTester::compareSolution(equationDiscretized4, referenceSolutionQuadratic, 1e-13);
  StiffnessMatrixTester::compareSolution(equationDiscretized5, referenceSolutionQuadratic, 1e-13);
  StiffnessMatrixTester::compareSolution(equationDiscretized6, referenceSolutionQuadratic, 1e-13);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect2DLinear)
{
  std::string pythonConfig = R"(
# Laplace 2D, Neumann BC

nx = 3
ny = nx

neumann_bc = []

# bottom, top
if True:
  for i in range(nx):
    neumann_bc += [
        {"element": i, "constantValue": -1.0, "face": "1-"},               # bottom
        {"element": (ny-1)*nx + i, "constantValue": 1.0, "face": "1+"},   # top
      ]

# left, right
if False:
  for j in range(ny):
    neumann_bc += [
        {"element": j*nx, "constantValue": -1, "face": "0-"},            # left
        {"element": j*nx + (nx-1), "constantValue": 1, "face": "0+"},   # top
      ]

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "nodePositions": [[float(i)/nx, float(j)/ny] for j in range(ny+1) for i in range(nx+1)],
    "elements": [[j*(nx+1)+i, j*(nx+1)+i+1, (j+1)*(nx+1)+i, (j+1)*(nx+1)+i+1] for j in range(ny) for i in range(nx)],
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : []
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized1(settings);
  equationDiscretized1.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized2(settings);
  equationDiscretized2.run();

  // unstructured
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized3(settings);
  equationDiscretized3.run();

  std::vector<double> referenceSolutionLinear = {
      0, 0, 0, 0, 1./3, 1./3, 1./3, 1./3, 2./3, 2./3, 2./3, 2./3, 1, 1, 1, 1
  };

  // for unstructured grids the order of the dofs is different
  std::vector<double> referenceSolutionLinearUnstructured = {
      0, 0, 1./3, 1./3, 0, 1./3, 0, 1./3, 2./3, 2./3, 2./3, 2./3, 1,  1,  1,  1
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionLinear);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionLinear);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionLinearUnstructured);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect2DQuadratic)
{
  std::string pythonConfig = R"(
# Laplace 2D, Neumann BC

nx = 3
ny = nx

neumann_bc = []

# bottom, top
if True:
  for i in range(nx):
    neumann_bc += [
        {"element": i, "constantValue": -1.0, "face": "1-"},               # bottom
        {"element": (ny-1)*nx + i, "constantValue": 1.0, "face": "1+"},   # top
      ]

# left, right
if False:
  for j in range(ny):
    neumann_bc += [
        {"element": j*nx, "constantValue": -1, "face": "0-"},            # left
        {"element": j*nx + (nx-1), "constantValue": 1, "face": "0+"},   # top
      ]

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "nodePositions": [[float(i)/(2*nx), float(j)/(2*ny)] for j in range(2*ny+1) for i in range(2*nx+1)],
    "elements": [[2*j*(2*nx+1)+2*i, 2*j*(2*nx+1)+2*i+1, 2*j*(2*nx+1)+2*i+2, (2*j+1)*(2*nx+1)+2*i, (2*j+1)*(2*nx+1)+2*i+1, (2*j+1)*(2*nx+1)+2*i+2, (2*j+2)*(2*nx+1)+2*i, (2*j+2)*(2*nx+1)+2*i+1, (2*j+2)*(2*nx+1)+2*i+2] for j in range(ny) for i in range(nx)],
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : []
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized1(settings);
  equationDiscretized1.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized2(settings);
  equationDiscretized2.run();

  // unstructured
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settings);
  equationDiscretized3.run();

  std::vector<double> referenceSolutionQuadratic = {
       0,  0, 0, 0, 0, 0, 0, 1./6,  1./6,  1./6,  1./6,  1./6,  1./6,  1./6,  1./3,  1./3,  1./3,  1./3,  1./3,  1./3,  1./3, 0.5,0.5,0.5,0.5,0.5,0.5,0.5, 2./3,  2./3,  2./3,  2./3,  2./3,  2./3,  2./3,  5./6,  5./6,  5./6,  5./6,  5./6,  5./6,  5./6,  1,  1,  1,  1,  1,  1,  1
  };

  // for unstructured grids the order of the dofs is different
  std::vector<double> referenceSolutionUnstructured = {
      0, 0, 0,  1./6,  1./6,  1./6,  1./3,  1./3,  1./3,  0, 0,  1./6,  1./6,  1./3,  1./3,  0, 0,  1./6,  1./6,  1./3,  1./3,  0.5,  0.5,  0.5,  2./3,  2./3,  2./3,  0.5,  0.5,  2./3,  2./3,  0.5,  0.5,  2./3,  2./3,  5./6,  5./6,  5./6,  1,  1,  1,  5./6,  5./6,  1,  1,  5./6,  5./6,  1,  1
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionQuadratic, 1e-13);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionQuadratic, 1e-13);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionUnstructured, 1e-13);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect2DHermite)
{
  std::string pythonConfig = R"(
# Laplace 2D, Neumann BC

nx = 3
ny = nx

neumann_bc = []

# bottom, top
if True:
  for i in range(nx):
    neumann_bc += [
        {"element": i, "constantValue": -1.0, "face": "1-"},               # bottom
        {"element": (ny-1)*nx + i, "constantValue": 1.0, "face": "1+"},   # top
      ]

# left, right
if False:
  for j in range(ny):
    neumann_bc += [
        {"element": j*nx, "constantValue": -1, "face": "0-"},            # left
        {"element": j*nx + (nx-1), "constantValue": 1, "face": "0+"},   # top
      ]

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "nodePositions": [[float(i)/nx, float(j)/ny] for j in range(ny+1) for i in range(nx+1)],
    "setHermiteDerivatives": True,
    "elements": [[j*(nx+1)+i, j*(nx+1)+i+1, (j+1)*(nx+1)+i, (j+1)*(nx+1)+i+1] for j in range(ny) for i in range(nx)],
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": neumann_bc,        # Neumann BC are always interpolated using Lagrange ansatz functions with one dof per node (not Hermite), even if the solution uses Hermite ansatz functions
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized1(settings);
  equationDiscretized1.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized2(settings);
  equationDiscretized2.run();

  // unstructured
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<2>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settings);
  equationDiscretized3.run();

  std::vector<double> referenceSolutionHermite = {
       0, 0, 1./3, 0, 0, 0, 1./3, 0, 0, 0, 1./3, 0, 0, 0, 1./3, 0, 1./3, 0, 1./3, 0, 1./3, 0, 1./3, 0, 1./3, 0, 1./3, 0, 1./3, 0, 1./3, 0, 2./3, 0, 1./3, 0, 2./3, 0, 1./3, 0, 2./3, 0, 1./3, 0, 2./3, 0, 1./3, 0, 1, 0, 1./3, 0, 1, 0, 1./3, 0, 1, 0, 1./3, 0, 1, 0, 1./3, 0
  };

  // for unstructured grids the order of the dofs is different
  std::vector<double> referenceSolutionUnstructured = {
      0, 0.406837, 0.281594, -2.24366, 0.030306, 0.081652, 0.0835685, -0.450449, 0.378492, 0.0185906, 0.0594175, -0.450268, 0.379885, 0.00320504, 0.0153181, -0.0873346, 0.030306, -0.081652, 0.0835685, 0.450449, 0.379885, -0.00320504, 0.0153181, 0.0873346, 8.52526e-15, -0.406837, 0.281594, 2.24366, 0.378492, -0.0185906, 0.0594175, 0.450268, 0.716161, -0.0185906, 0.0594175, -0.450268, 0.714768, -0.00320504, 0.0153181, -0.0873346, 0.714768, 0.00320504, 0.0153181, 0.0873346, 0.716161, 0.0185906, 0.0594175, 0.450268, 1.09465, -0.406837, 0.281594, -2.24366, 1.06435, -0.081652, 0.0835685, -0.450449, 1.06435, 0.081652, 0.0835685, 0.450449, 1.09465, 0.406837, 0.281594, 2.24366
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionHermite, 1e-12);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionHermite, 1e-12);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionUnstructured, 1e-5);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect3DLinearQuadratic)
{

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
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace_structured", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}

)";

   std::string pythonConfigUnstructured = R"(
# Laplace 3D, Neumann BC
import sys

# 3D laplace problem

nx = 2
ny = 2
nz = 5

# Neumann boundary conditions
bc = []

# global boundary conditions
for j in range(int(ny)):
  for i in range(int(nx)):
    element_no = int(j*nx + i)

    # z- plane (bottom)
    bc_value = -1.0   # negative = inflow

    bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

    # z+ plane (top)
    element_no += int((nz-1)*nx*ny)
    bc_value = 1.0
    bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

n_elements = [nx,ny,nz]

# node positions
node_positions = []
for k in range(int(nz+1)):
  for j in range(int(ny+1)):
    for i in range(int(nx+1)):
      element_no = int(j*nx + i)
      node_positions.append([i * 1.0/nx, j * 1.0/ny, k * 1.0/nz])

# elements
elements = []
for k in range(int(nz)):
  for j in range(int(ny)):
    for i in range(int(nx)):
      elements.append([k*(nx+1)*(ny+1) + j*(nx+1) + i, k*(nx+1)*(ny+1) + j*(nx+1) + i+1, k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i, k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1,
                       (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i, (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i+1, (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i, (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1,
      ])

config = {
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "inputMeshIsGlobal": True,
    "nodePositions": node_positions,
    "elements": elements,
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace_linear_unstructured", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}

)";

  DihuContext settings(argc, argv, pythonConfig);

  // regular fixed
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized1(settings);
  equationDiscretized1.run();

  // structured deformable
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized2(settings);
  equationDiscretized2.run();

  // unstructured
  DihuContext settingsUnstructured(argc, argv, pythonConfigUnstructured);
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settingsUnstructured);
  equationDiscretized3.run();

  std::vector<double> referenceSolutionLinear = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5
  };

  // for unstructured grids the order of the dofs is different
  std::vector<double> referenceSolutionLinearUnstructured = {
      0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0, 0, 0.2, 0.2, 0, 0, 0.2, 0.2, 0, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionLinear, 1e-5);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionLinear, 1e-5);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionLinearUnstructured, 1e-5);

  // quadratic

  // regular fixed
  LOG(DEBUG) << "======================== quadratic regular fixed =======================";
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized4(settings);
  equationDiscretized4.run();

  // structured deformable
  LOG(DEBUG) << "======================== quadratic structured deformable =======================";
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized5(settings);
  equationDiscretized5.run();

  std::vector<double> referenceSolutionQuadratic = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
  };

  //StiffnessMatrixTester::compareSolution(equationDiscretized4, referenceSolutionQuadratic, 1e-5);
  StiffnessMatrixTester::compareSolution(equationDiscretized5, referenceSolutionQuadratic, 1e-5);

  // Hermite

  // regular fixed
  LOG(DEBUG) << "======================== Hermite regular fixed =======================";
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized6(settings);
  equationDiscretized6.run();

  // structured deformable
  LOG(DEBUG) << "======================== Hermite structured deformable =======================";
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized7(settings);
  equationDiscretized7.run();

  std::vector<double> referenceSolutionHermite = {
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0, 5, 0, 0, 0, 1, 0, 0, 0
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized6, referenceSolutionHermite, 1e-5);
  StiffnessMatrixTester::compareSolution(equationDiscretized7, referenceSolutionHermite, 1e-5);

}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect3DQuadraticUnstructured)
{
   std::string pythonConfigUnstructured = R"(
# Laplace 3D, Neumann BC
import sys

# 3D laplace problem

nx = 2
ny = 2
nz = 5

# Neumann boundary conditions
bc = []

# global boundary conditions
for j in range(int(ny)):
  for i in range(int(nx)):
    element_no = int(j*nx + i)

    # z- plane (bottom)
    bc_value = -1.0   # negative = inflow

    bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

    # z+ plane (top)
    element_no += int((nz-1)*nx*ny)
    bc_value = 1.0
    bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

n_elements = [nx,ny,nz]

# node positions
node_positions = []
for k in range(int(2*nz+1)):
  for j in range(int(2*ny+1)):
    for i in range(int(2*nx+1)):
      node_positions.append([i * 1.0/(2*nx), j * 1.0/(2*ny), k * 1.0/(2*nz)])

# elements
elements = []
for k in range(int(nz)):
  for j in range(int(ny)):
    for i in range(int(nx)):
      elements.append([2*k*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i, 2*k*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i+1,  2*k*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i+2,
                       2*k*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i, 2*k*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i+1,  2*k*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i+2,
                       2*k*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i, 2*k*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i+1,  2*k*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i+2,
                       (2*k+1)*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i, (2*k+1)*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i+1,  (2*k+1)*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i+2,
                       (2*k+1)*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i, (2*k+1)*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i+1,  (2*k+1)*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i+2,
                       (2*k+1)*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i, (2*k+1)*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i+1,  (2*k+1)*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i+2,
                       (2*k+2)*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i, (2*k+2)*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i+1,  (2*k+2)*(2*nx+1)*(2*ny+1) + 2*j*(2*nx+1) + 2*i+2,
                       (2*k+2)*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i, (2*k+2)*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i+1,  (2*k+2)*(2*nx+1)*(2*ny+1) + (2*j+1)*(2*nx+1) + 2*i+2,
                       (2*k+2)*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i, (2*k+2)*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i+1,  (2*k+2)*(2*nx+1)*(2*ny+1) + (2*j+2)*(2*nx+1) + 2*i+2
                       ])

config = {
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "inputMeshIsGlobal": True,
    "nodePositions": node_positions,
    "elements": elements,
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace_quadratic_unstructured", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}

)";

  // unstructured
  DihuContext settingsUnstructured(argc, argv, pythonConfigUnstructured);
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settingsUnstructured);
  equationDiscretized3.run();

  // for unstructured grids the order of the dofs is different
  std::vector<double> referenceSolutionLinearUnstructured = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.8, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1, 1, 1, 1, 1, 1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1, 1, 1, 1, 1, 1, 0.9, 0.9, 0.9, 0.9, 1, 1, 1, 1
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionLinearUnstructured, 1e-5);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect3DHermiteUnstructured)
{
   std::string pythonConfigUnstructured = R"(
# Laplace 3D, Neumann BC
import sys

# 3D laplace problem

nx = 2
ny = 2
nz = 5

# Neumann boundary conditions
bc = []

# global boundary conditions
for j in range(int(ny)):
  for i in range(int(nx)):
    element_no = int(j*nx + i)

    # z- plane (bottom)
    bc_value = -1.0   # negative = inflow

    bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})

    # z+ plane (top)
    element_no += int((nz-1)*nx*ny)
    bc_value = 1.0
    bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

n_elements = [nx,ny,nz]

# node positions
node_positions = []
for k in range(int(nz+1)):
  for j in range(int(ny+1)):
    for i in range(int(nx+1)):
      node_positions.append([i * 1.0/nx, j * 1.0/ny, k * 1.0/nz])

# elements
elements = []
for k in range(int(nz)):
  for j in range(int(ny)):
    for i in range(int(nx)):

      elements.append([k*(nx+1)*(ny+1) + j*(nx+1) + i, k*(nx+1)*(ny+1) + j*(nx+1) + i+1, k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i, k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1, \
                            (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i, (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i+1, (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i, (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1])

config = {
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "inputMeshIsGlobal": True,
    "nodePositions": node_positions,
    "elements": elements,
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "sor",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace_hermite_unstructured", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}

)";

  // unstructured
  DihuContext settingsUnstructured(argc, argv, pythonConfigUnstructured);
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<3>,
    BasisFunction::Hermite,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized(settingsUnstructured);
  equationDiscretized.run();

  // for unstructured grids the order of the dofs is different
  std::vector<double> referenceSolutionHermiteUnstructured = {
      0, -0.495451, -0.49651, 6.43463, -0.349013, 3.88378, 3.89475, -37.0353, -0.179265, -0.435614, 0.69288, 3.16961, 0.344894, 2.42984, -1.85321, -17.0209, -0.179347, 0.693795, -0.43699, 3.17497, 0.34577, -1.86114, 2.44285, -17.0907, -0.138322, -0.0405356, -0.0408618, 1.00509, 0.217959, 0.189338, 0.193204, -5.3168, 0.16836, -0.071325, -0.071753, 0.552223, -0.101576, 0.812085, 0.814034, -6.18138, 0.152111, -0.0288353, -0.0172367, 0.219415, -0.0244468, 0.34471, 0.149262, -2.62629, 0.152125, -0.0169803, -0.0292689, 0.22138, -0.0241534, 0.146859, 0.346633, -2.63798, 0.147175, -0.0033467, -0.00348111, 0.0727128, -0.0120092, 0.024662, 0.0256294, -0.970351, -0.328558, -1.73757, 1.51358, 8.71324, 1.01352, 7.02171, -5.70884, -39.4292, -0.201597, -0.866968, 0.0977598, 1.22472, 0.451829, 2.70815, -0.554999, -6.76122, 0.149705, -0.0343287, 0.01695, 0.308237, 0.025364, 0.454246, -0.267125, -3.96883, 0.150075, 0.00833965, 0.00709421, 0.0906284, -0.0146, -0.0535205, -0.113083, -1.35005, -0.329031, 1.51735, -1.74113, 8.73652, 1.01747, -5.74012, 7.05344, -39.6511, -0.201727, 0.0991285, -0.867562, 1.23153, 0.452865, -0.563675, 2.71706, -6.81031, 0.149595, 0.0183041, -0.0351215, 0.316313, 0.0260524, -0.273627, 0.459457, -4.01445, 0.150009, 0.00749827, 0.00787043, 0.0919097, -0.0140988, -0.113498, -0.0495052, -1.35116, -0.277546, -1.08893, -1.08993, -3.37205, 0.701887, 3.21992, 3.22475, 9.13765, 0.154654, 0.0201825, 0.0198853, 0.119052, -0.0341112, -0.218288, -0.218912, -1.68434, 0.354081, -0.0247758, -0.0247433, 0.219604, -0.0355886, 0.307951, 0.30833, -2.46117, 0.34954, -0.0128317, -0.0025973, 0.103192, -0.00777113, 0.138016, 0.0468795, -1.12709, 0.349531, -0.00237932, -0.0128697, 0.104147, -0.00765516, 0.0461305, 0.137598, -1.12295, 0.348104, -0.000664865, -0.000791976, 0.0483819, -0.00445587, 0.00986536, 0.00964765, -0.477708, 0.347122, -0.0170531, 0.0123475, 0.160376, 0.0132641, 0.192239, -0.129849, -1.75699, 0.348158, 0.00147705, 0.00508855, 0.0684468, -0.00357497, -0.00918824, -0.0591085, -0.697664, 0.347178, 0.0118969, -0.0170393, 0.160634, 0.0140629, -0.132765, 0.195311, -1.76988, 0.348134, 0.00487778, 0.00117169, 0.0682855, -0.00348443, -0.0583755, -0.00867446, -0.694652, 0.349243, 0.00833332, 0.00869989, 0.0922472, -0.0121117, -0.0982075, -0.0991534, -0.915439, 0.549254, -0.00988294, -0.00970342, 0.11888, -0.019973, 0.150793, 0.150703, -1.08045, 0.548152, -0.00609839, 0.00108792, 0.0594731, -0.00668321, 0.0618669, 0.0390271, -0.527423, 0.548138, 0.00114646, -0.00600757, 0.0590501, -0.00647558, 0.0382294, 0.0609433, -0.515241, 0.548263, -0.000575121, -0.000602484, 0.033749, -0.0043158, 0.00766694, 0.00690726, -0.229796, 0.546042, -0.00967511, 0.00916097, 0.0864894, 0.00264034, 0.0743068, -0.0474516, -0.838166, 0.547265, -0.00187614, 0.00399411, 0.0507699, -0.00401416, -0.0128517, -0.0265187, -0.323011, 0.546042, 0.00935924, -0.00966486, 0.0872386, 0.00352116, -0.0521058, 0.0778085, -0.84587, 0.547253, 0.00395666, -0.00199061, 0.0506206, -0.00415121, -0.0272813, -0.0142736, -0.3199, 0.547239, 0.00576412, 0.00566091, 0.0844041, -0.00810736, -0.0534196, -0.0519008, -0.365263, 0.743491, -0.00108863, -0.000937559, 0.1298, -0.0283731, 0.112339, 0.113811, -0.105561, 0.74576, -0.00653243, 0.0110359, 0.0783566, -0.0188302, 0.0158994, 0.100156, -0.0621963, 0.745715, 0.0111441, -0.00646536, 0.0774174, -0.0188016, 0.098559, 0.0145459, -0.0508119, 0.749358, -0.00059729, -0.000423098, 0.0498373, -0.011661, 0.00514494, 0.00619304, 0.0169457, 0.740782, -0.0197251, 0.0202853, 0.118573, -0.0209067, -0.050852, 0.0713799, -0.166978, 0.745076, -0.0120117, 0.0055416, 0.0730212, -0.0172976, -0.0824509, 0.00827695, 0.0719195, 0.740654, 0.0207217, -0.0204753, 0.121141, -0.0204413, 0.0696096, -0.0491709, -0.172743, 0.745117, 0.0056067, -0.0116294, 0.0730088, -0.0175745, 0.0076752, -0.0831276, 0.0664507, 0.742112, -0.00140029, -0.00178678, 0.107925, -0.0201842, -0.0522637, -0.0508324, 0.31155, 1.18658, -1.24234, -1.24339, 5.23629, 0.768471, -4.04349, -4.04583, 19.3078, 1.09447, -0.0129252, -0.82185, 0.175796, 0.43189, -0.059938, -2.45772, 0.903187, 1.09441, -0.821513, -0.0146078, 0.188585, 0.431741, -2.45889, -0.0675975, 0.963794, 1.03673, 0.000874937, 0.000534494, 0.111551, 0.224281, 0.00735688, 0.00863065, 0.614291, 1.1822, 1.19948, -1.19494, -4.66283, 0.752967, 3.84913, -3.80867, -16.3803, 1.09388, 0.823809, 0.0147377, 0.169142, 0.432848, 2.47826, 0.0874204, 0.982568, 1.18179, -1.19278, 1.19589, -4.63529, 0.751656, -3.80142, 3.83527, -16.2816, 1.09396, 0.0152444, 0.823906, 0.175166, 0.433212, 0.0886698, 2.48136, 0.999248, 1.18557, 1.242, 1.24165, 5.19698, 0.774308, 4.08317, 4.08312, 19.422
  };

  StiffnessMatrixTester::compareSolution(equationDiscretized, referenceSolutionHermiteUnstructured, 5e-5);
}

}  // namespace

