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

TEST(LaplaceTest, MatrixIsCorrect1DSmall)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 5
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "FiniteElementMethod": {
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();

  std::vector<double> referenceMatrix = {
    1, 0, 0, 0, 0, 0,
    0, -2.5, 1.25, 0, 0, 0,
    0, 1.25, -2.5, 1.25, 0, 0,
    0, 0, 1.25, -2.5, 1.25, 0, 
    0, 0, 0, 1.25, -2.5, 0,
    0, 0, 0, 0, 0, 1
  };
  std::vector<double> referenceRhs = {
    1, -1.25, 0, 0, 0, 0
  };
  std::map<int, double> dirichletBC = {{0, 1.0}, {5,0.0}};
  
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized, dirichletBC);
}

TEST(LaplaceTest, MatrixIsCorrect1DBig)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 10
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 2.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 0.5,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();

  std::vector<double> referenceMatrix = {
    1,    0,    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, -40, 20, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 20, -40, 20, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 20, -40, 20, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 20, -40, 20, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 20, -40, 20, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 20, -40, 20, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 20, -40, 20, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 20, -40, 20, 0,
    0, 0, 0, 0, 0, 0, 0, 0,    20, -40, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  };
  std::vector<double> referenceRhs = {
    1, -20, 0, 0, 0, 0, 0, 0, 0, -40, 2
  };
  std::map<int, double> dirichletBC = {{0, 1.0}, {10,2.0}};
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized, dirichletBC);
}

TEST(LaplaceTest, MatrixIsCorrect1DStencils)
{
  std::string pythonConfig = R"(
# Laplace 1D
n = 10
    
# no boundary conditions
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": n,
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();

  std::vector<double> referenceMatrix(121, 0.0);
  std::array<double, 2> stencil = {-1.0, 1.0};
    // stencil for -Δu in 1D: [1 _-2_ 1] (element contribution: [_-1_ 1])
  
  // fill with stencil values
  auto matrixIndex = [](int x, int y){return y*11 + x;};
  
  // loop over elements
  for (int i=0; i<10; i++)
  {
    // nodes:
    // 0--1
    int node0 = i;
    int node1 = i+1;
    
    // add contribution from node 0 to all nodes
    referenceMatrix[matrixIndex(node0, node0)] += stencil[0];
    referenceMatrix[matrixIndex(node0, node1)] += stencil[1];
    
    // add contribution from node 1 to all nodes
    referenceMatrix[matrixIndex(node1, node0)] += stencil[1];
    referenceMatrix[matrixIndex(node1, node1)] += stencil[0];
  }
  
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}


}  // namespace

