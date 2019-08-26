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

}  // namespace

