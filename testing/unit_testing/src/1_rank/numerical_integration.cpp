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
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);

  equationDiscretized1.run();

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

TEST(NumericalIntegrationTest, StiffnessMatrixIsCorrect1DBCBackwards)
{
  std::string pythonConfig1 = R"(
# Laplace 1D
n = 5

# boundary conditions
bc = {}
bc[0] = 1.0
bc[-1] = 5.0

config = {
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
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);

  equationDiscretized1.run();

  std::vector<double> referenceMatrix = {
    1, 0, 0, 0, 0, 0,
    0, -2.5, 1.25, 0, 0, 0,
    0, 1.25, -2.5, 1.25, 0, 0,
    0, 0, 1.25, -2.5, 1.25, 0,
    0, 0, 0, 1.25, -2.5, 0,
    0, 0, 0, 0, 0, 1
  };
  std::vector<double> referenceRhs = {
    1, -1.25, 0, 0, -6.25, 5
  };
  std::map<int, double> dirichletBC = {{0, 1.0}, {5,5.0}};

  StiffnessMatrixTester::compareMatrix(equationDiscretized1, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized1, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized1, dirichletBC);
}
/*

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
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  equationDiscretized1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // integration order: lagrange basis order 1 => ∇phi_i*∇phi_j order 2 => 2 gauss points => p_exact=3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  equationDiscretized2.run();
  
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
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  equationDiscretized1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // integration order: lagrange basis order 1 => ∇phi_i*∇phi_j order 2 => 2 gauss points => p_exact=3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  equationDiscretized2.run();
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized2);
}*/
/*
 * // the following tests are commented out because the take very long to compile, they should work, however
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
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  equationDiscretized1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // Gauss integration order 3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settings2);
  
  equationDiscretized3.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized3);
  
  // Gauss integration order 4
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<4>,
    Equation::Static::Laplace
  > equationDiscretized4(settings2);
  
  equationDiscretized4.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized4);
  
  // Gauss integration order 5
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<5>,
    Equation::Static::Laplace
  > equationDiscretized5(settings2);
  
  equationDiscretized5.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized5);
  
  // Gauss integration order 6
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<6>,
    Equation::Static::Laplace
  > equationDiscretized6(settings2);
  
  equationDiscretized6.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized6);
  
  // Gauss integration order 7
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<7>,
    Equation::Static::Laplace
  > equationDiscretized7(settings2);
  
  equationDiscretized7.run(); 
  
  // Gauss integration order 8
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<8>,
    Equation::Static::Laplace
  > equationDiscretized8(settings2);
  
  equationDiscretized8.run(); 
  
  // Gauss integration order 10
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<10>,
    Equation::Static::Laplace
  > equationDiscretized10(settings2);
  
  equationDiscretized10.run(); 
  
  // Gauss integration order 12
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<12>,
    Equation::Static::Laplace
  > equationDiscretized12(settings2);
  
  equationDiscretized12.run(); 
  

  // Gauss integration order 16
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<16>,
    Equation::Static::Laplace
  > equationDiscretized16(settings2);
  
  equationDiscretized16.run(); 
  

  // Gauss integration order 20
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<20>,
    Equation::Static::Laplace
  > equationDiscretized20(settings2);
  
  equationDiscretized20.run(); 
  

  // Gauss integration order 24
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<24>,
    Equation::Static::Laplace
  > equationDiscretized24(settings2);
  
  equationDiscretized24.run(); 
  
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized7);
  
  // Gauss integration order 64
//  FiniteElementMethod<
//    Mesh::StructuredDeformableOfDimension<3>,
//    BasisFunction::LagrangeOfOrder<>,
//    Quadrature::Gauss<64>,
//    Equation::Static::Laplace
  //> equationDiscretized64(settings2);
  
  //Computation computation64(settings2, equationDiscretized64);
  //computation64.run(); 
  
  //StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized64);
  
}

TEST(NumericalIntegrationTest, NewtonCotesIntegrationHigherOrderWorks)
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
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  equationDiscretized1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // NewtonCotes integration order 3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::NewtonCotes<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settings2);
  
  equationDiscretized3.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized3);
  
  // NewtonCotes integration order 4
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::NewtonCotes<4>,
    Equation::Static::Laplace
  > equationDiscretized4(settings2);
  
  equationDiscretized4.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized4);
  
  // NewtonCotes integration order 5
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::NewtonCotes<5>,
    Equation::Static::Laplace
  > equationDiscretized5(settings2);
  
  equationDiscretized5.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized5);
  
  // NewtonCotes integration order 6
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::NewtonCotes<6>,
    Equation::Static::Laplace
  > equationDiscretized6(settings2);
  
  equationDiscretized6.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized6);
  
  // NewtonCotes integration order 7
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::NewtonCotes<7>,
    Equation::Static::Laplace
  > equationDiscretized7(settings2);
  
  equationDiscretized7.run(); 
  
  // NewtonCotes integration order 8
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::NewtonCotes<8>,
    Equation::Static::Laplace
  > equationDiscretized8(settings2);
  
  equationDiscretized8.run(); 
}

TEST(NumericalIntegrationTest, ClenshawCurtisIntegrationHigherOrderWorks)
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
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized1(settings1);
  
  equationDiscretized1.run();

  DihuContext settings2(argc, argv, pythonConfig);
  
  // ClenshawCurtis integration order 3
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::ClenshawCurtis<3>,
    Equation::Static::Laplace
  > equationDiscretized3(settings2);
  
  equationDiscretized3.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized3);
  
  // ClenshawCurtis integration order 4
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::ClenshawCurtis<4>,
    Equation::Static::Laplace
  > equationDiscretized4(settings2);
  
  equationDiscretized4.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized4);
  
  // ClenshawCurtis integration order 5
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::ClenshawCurtis<5>,
    Equation::Static::Laplace
  > equationDiscretized5(settings2);
  
  equationDiscretized5.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized5);
  
  // ClenshawCurtis integration order 6
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::ClenshawCurtis<6>,
    Equation::Static::Laplace
  > equationDiscretized6(settings2);
  
  equationDiscretized6.run(); 
  
  StiffnessMatrixTester::checkEqual(equationDiscretized1, equationDiscretized6);
  
  // ClenshawCurtis integration order 7
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::ClenshawCurtis<7>,
    Equation::Static::Laplace
  > equationDiscretized7(settings2);
  
  equationDiscretized7.run(); 
}
*/
};


