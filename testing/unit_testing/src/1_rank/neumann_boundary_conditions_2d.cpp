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

}  // namespace

