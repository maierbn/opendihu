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

/*
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
    "preconditionerType": "none",
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
    "preconditionerType": "none",
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

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionQuadratic);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionQuadratic);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionUnstructured);
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
    "preconditionerType": "none",
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

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionHermite, 1e-13);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionHermite, 1e-12);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionUnstructured, 1e-5);
}

*/
/*
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
    "preconditionerType": "none",
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
    "preconditionerType": "none",
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
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized4(settings);
  equationDiscretized4.run();

  // structured deformable
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


  StiffnessMatrixTester::compareSolution(equationDiscretized4, referenceSolutionQuadratic, 1e-5);
  StiffnessMatrixTester::compareSolution(equationDiscretized5, referenceSolutionQuadratic, 1e-5);
}*/


TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect3DQuadratic)
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
    "preconditionerType": "none",
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
/*
TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect3DQuadratic)
{
  std::string pythonConfig = R"(
# Laplace 3D, Neumann BC

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
    "nodePositions": [[float(i)/(2*nx), float(j)/(2*ny)] for j in range(2*ny+1) for i in range(2*nx+1)],    "elements": [[2*j*(2*nx+1)+2*i, 2*j*(2*nx+1)+2*i+1, 2*j*(2*nx+1)+2*i+2, (2*j+1)*(2*nx+1)+2*i, (2*j+1)*(2*nx+1)+2*i+1, (2*j+1)*(2*nx+1)+2*i+2, (2*j+2)*(2*nx+1)+2*i, (2*j+2)*(2*nx+1)+2*i+1, (2*j+2)*(2*nx+1)+2*i+2] for j in range(ny) for i in range(nx)],
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0},
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
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

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionQuadratic);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionQuadratic);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionUnstructured);
}

TEST(NeumannBoundaryConditionTest, NeumannBoundaryConditionIsCorrect3DHermite)
{
  std::string pythonConfig = R"(
# Laplace 3D, Neumann BC

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
    "preconditionerType": "none",
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

  StiffnessMatrixTester::compareSolution(equationDiscretized1, referenceSolutionHermite, 1e-13);
  StiffnessMatrixTester::compareSolution(equationDiscretized2, referenceSolutionHermite, 1e-12);
  StiffnessMatrixTester::compareSolution(equationDiscretized3, referenceSolutionUnstructured, 1e-5);
}

*/
}  // namespace

