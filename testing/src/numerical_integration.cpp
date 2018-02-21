#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "utility/petsc_utility.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"

namespace SpatialDiscretization
{
  
TEST(NumericalIntegrationTest, StiffnessMatrixIsCorrect1D)
{
  std::string pythonConfig1 = R"(
# Laplace 1D
n = 5
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  }
}
)";

  DihuContext settings1(argc, argv, pythonConfig1);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  Computation computation(settings1, equationDiscretized1);
  computation.run();

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
  
  StiffnessMatrixTester::compareMatrix(equationDiscretized1, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized1, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized1, dirichletBC);
}


TEST(NumericalIntegrationTest, StiffnessMatrixIsCorrect2D)
{
  std::string pythonConfig = R"(
# Laplace 2D
n = 2
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "prefactor": 6,
    "nElements": [n,n],
    "physicalExtent": [4.0,4.0],
    #"DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  }
}
)";

  DihuContext settings1(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  Computation computation1(settings1, equationDiscretized1);
  computation1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // integration order: lagrange basis order 1 => ∇phi_i*∇phi_j order 2 => 2 gauss points => p_exact=3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  Computation computation2(settings2, equationDiscretized2);
  computation2.run();
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized2);
}


TEST(NumericalIntegrationTest, StiffnessMatrixIsCorrect3D)
{
  std::string pythonConfig = R"(
# Laplace 3D
n = 2
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "prefactor": 12,
    "nElements": [n,n,n],
    "physicalExtent": [3.0,3.0,3.0],
    #"DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  }
}
)";

  DihuContext settings1(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  Computation computation1(settings1, equationDiscretized1);
  computation1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // integration order: lagrange basis order 1 => ∇phi_i*∇phi_j order 2 => 2 gauss points => p_exact=3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  Computation computation2(settings2, equationDiscretized2);
  computation2.run();
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized2);
}

TEST(NumericalIntegrationTest, GaussIntegrationHigherOrderWorks)
{
  std::string pythonConfig = R"(
# Laplace 3D
n = 2
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "prefactor": 12,
    "nElements": [n,n,n],
    "physicalExtent": [2.0,2.0,2.0],
    #"DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
  }
}
)";

  DihuContext settings1(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  Computation computation1(settings1, equationDiscretized1);
  computation1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // Gauss integration order 3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settings2);
  
  Computation computation3(settings2, equationDiscretized3);
  computation3.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized3);
  
  // Gauss integration order 4
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<4>,
    Equation::Static::Laplace
  > equationDiscretized4(settings2);
  
  Computation computation4(settings2, equationDiscretized4);
  computation4.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized4);
  
  // Gauss integration order 5
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<5>,
    Equation::Static::Laplace
  > equationDiscretized5(settings2);
  
  Computation computation5(settings2, equationDiscretized5);
  computation5.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized5);
  
  // Gauss integration order 6
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<6>,
    Equation::Static::Laplace
  > equationDiscretized6(settings2);
  
  Computation computation6(settings2, equationDiscretized6);
  computation6.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized6);
  
  // Gauss integration order 7
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Integrator::Gauss<7>,
    Equation::Static::Laplace
  > equationDiscretized7(settings2);
  
  Computation computation7(settings2, equationDiscretized7);
  computation7.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized7);
  
}

};

