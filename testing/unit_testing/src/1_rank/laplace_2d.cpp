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

TEST(LaplaceTest, MatrixIsCorrect2DSmall)
{
  std::string pythonConfig = R"(
# Laplace 2D,  3x4 nodes
n = 2.
m = 3.

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = x
  i2 = (n+1)*m + i
  bc[int(i2)] = x
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

#bc = {0: 0, 1: 1./3, 2: 2./3,
#      9: 0, 10: 1./3, 11: 2./3}

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "physicalExtent": [6.0, 9.0],
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();

  std::vector<double> referenceMatrix = {
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0,-12/9.,3/9., 0,1.5/9.,3/9., 0, 0, 0, 0,
    0, 0, 0,3/9.,-24/9.,3/9.,3/9.,3/9., 3/9., 0, 0, 0,
    0, 0, 0, 0,3/9.,-12/9.,0,3/9.,1.5/9., 0, 0, 0,
    0, 0, 0,1.5/9.,3/9.,0,-12/9.,3/9., 0, 0, 0, 0,
    0, 0, 0,3/9.,3/9.,3/9.,3/9.,-24/9., 3/9., 0, 0, 0,
    0, 0, 0,0,3/9.,1.5/9.,0,3/9.,-12/9., 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  };
  std::vector<double> referenceRhs = {
    0, 1./3, 2./3, -1./9, -3./9, -2./9, -1./9, -3/9., -2./9, 0, 1./3, 2./3
  };
  std::map<int, double> dirichletBC = {{0, 0.0}, {1,1./3}, {2,2./3}, {9,0}, {10,1./3}, {11,2./3}};
  
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized, dirichletBC);
}

TEST(LaplaceTest, MatrixIsCorrect2DStencils)
{
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 4],
    "physicalExtent": [4.0, 4.0],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();

  std::vector<double> referenceMatrix(625, 0.0);
  
  // stencil for -Δu in 2D:     [1  1   1] (element contribution: [  1/6  1/3])
  //                        1/3*[1 _-8_ 1]                        [_-2/3_ 1/6]
  //                            [1  1   1]
  double stencil[2][2] = {{-2./3, 1./6}, {1./6, 1./3}};
  
  // fill with stencil values
  auto matrixIndex = [](int x, int y){return y*25 + x;};
  
  // loop over elements
  for (int x=0; x<4; x++)
  {
    for (int y=0; y<4; y++)
    {
      // nodes:
      // 2--3
      // 0--1
      int node0 = y*5 + x;
      int node1 = y*5 + x+1;
      int node2 = (y+1)*5 + x;
      int node3 = (y+1)*5 + x+1;
      
      // add contribution from node 0 to all nodes
      referenceMatrix[matrixIndex(node0, node0)] += stencil[0][0];
      referenceMatrix[matrixIndex(node0, node1)] += stencil[0][1];
      referenceMatrix[matrixIndex(node0, node2)] += stencil[1][0];
      referenceMatrix[matrixIndex(node0, node3)] += stencil[1][1];
      
      // add contribution from node 1 to all nodes
      referenceMatrix[matrixIndex(node1, node0)] += stencil[0][1];
      referenceMatrix[matrixIndex(node1, node1)] += stencil[0][0];
      referenceMatrix[matrixIndex(node1, node2)] += stencil[1][1];
      referenceMatrix[matrixIndex(node1, node3)] += stencil[1][0];
      
      // add contribution from node 2 to all nodes
      referenceMatrix[matrixIndex(node2, node0)] += stencil[1][0];
      referenceMatrix[matrixIndex(node2, node1)] += stencil[1][1];
      referenceMatrix[matrixIndex(node2, node2)] += stencil[0][0];
      referenceMatrix[matrixIndex(node2, node3)] += stencil[0][1];
      
      // add contribution from node 3 to all nodes
      referenceMatrix[matrixIndex(node3, node0)] += stencil[1][1];
      referenceMatrix[matrixIndex(node3, node1)] += stencil[1][0];
      referenceMatrix[matrixIndex(node3, node2)] += stencil[0][1];
      referenceMatrix[matrixIndex(node3, node3)] += stencil[0][0];
    }
  }
    
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}

TEST(LaplaceTest, StructuredDeformableMatrixIsCorrect2DStencils)
{
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 4],
    "physicalExtent": [4.0, 4.0],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  equationDiscretized.run();

  std::vector<double> referenceMatrix(625, 0.0);
  
  // stencil for -Δu in 2D:     [1  1   1] (element contribution: [  1/6  1/3])
  //                        1/3*[1 _-8_ 1]                        [_-2/3_ 1/6]
  //                            [1  1   1]
  double stencil[2][2] = {{-2./3, 1./6}, {1./6, 1./3}};
  
  // fill with stencil values
  auto matrixIndex = [](int x, int y){return y*25 + x;};
  
  // loop over elements
  for (int x=0; x<4; x++)
  {
    for (int y=0; y<4; y++)
    {
      // nodes:
      // 2--3
      // 0--1
      int node0 = y*5 + x;
      int node1 = y*5 + x+1;
      int node2 = (y+1)*5 + x;
      int node3 = (y+1)*5 + x+1;
      
      // add contribution from node 0 to all nodes
      referenceMatrix[matrixIndex(node0, node0)] += stencil[0][0];
      referenceMatrix[matrixIndex(node0, node1)] += stencil[0][1];
      referenceMatrix[matrixIndex(node0, node2)] += stencil[1][0];
      referenceMatrix[matrixIndex(node0, node3)] += stencil[1][1];
      
      // add contribution from node 1 to all nodes
      referenceMatrix[matrixIndex(node1, node0)] += stencil[0][1];
      referenceMatrix[matrixIndex(node1, node1)] += stencil[0][0];
      referenceMatrix[matrixIndex(node1, node2)] += stencil[1][1];
      referenceMatrix[matrixIndex(node1, node3)] += stencil[1][0];
      
      // add contribution from node 2 to all nodes
      referenceMatrix[matrixIndex(node2, node0)] += stencil[1][0];
      referenceMatrix[matrixIndex(node2, node1)] += stencil[1][1];
      referenceMatrix[matrixIndex(node2, node2)] += stencil[0][0];
      referenceMatrix[matrixIndex(node2, node3)] += stencil[0][1];
      
      // add contribution from node 3 to all nodes
      referenceMatrix[matrixIndex(node3, node0)] += stencil[1][1];
      referenceMatrix[matrixIndex(node3, node1)] += stencil[1][0];
      referenceMatrix[matrixIndex(node3, node2)] += stencil[0][1];
      referenceMatrix[matrixIndex(node3, node3)] += stencil[0][0];
    }
  }
    
  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}


}  // namespace

