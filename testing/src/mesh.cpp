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
  
TEST(MeshTest, ReadNodePositionsAreCorrect1D)
{
  // explicit mesh with node positions
  std::string pythonConfig = R"(
# Laplace 1D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "Meshes" : {
    "testMesh": {
      "nElements": [6],
      "nodeDimension": 1,
      "nodePositions": [0,1,2,3,4,5,6],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  std::vector<double> referenceNodePositions = {
    0,0,0,  1,0,0,  2,0,0,
    3,0,0,  4,0,0,  5,0,0,
    6,0,0,
  };
  
  Mesh::NodePositionsTester::compareNodePositions(settings, "testMesh", referenceNodePositions);
  
  // explicit mesh with automatically generated node positions
  std::string pythonConfig2 = R"(
# Laplace 1D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "Meshes" : {
    "testMesh": {
      "nElements": 6,
      "physicalExtent": 6,
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  Computation computation2(settings2, equationDiscretized2);
  computation2.run();
  
  Mesh::NodePositionsTester::compareNodePositions(settings2, "testMesh", referenceNodePositions);
}

TEST(MeshTest, ReadNodePositionsAreCorrect2D)
{
  // explicit mesh with node positions
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "Meshes" : {
    "testMesh": {
      "nElements": [2, 2],
      "nodeDimension": 2,
      "nodePositions": [0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  std::vector<double> referenceNodePositions = {
    0,0,0,  1,0,0,  2,0,0,
    0,1,0,  1,1,0,  2,1,0,
    0,2,0,  1,2,0,  2,2,0,
  };
  
  Mesh::NodePositionsTester::compareNodePositions(settings, "testMesh", referenceNodePositions);
  
  // explicit mesh with automatically generated node positions
  std::string pythonConfig2 = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "Meshes" : {
    "testMesh": {
      "nElements": [2, 2],
      "physicalExtent": [2, 2],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<1>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  Computation computation2(settings2, equationDiscretized2);
  computation2.run();
  
  Mesh::NodePositionsTester::compareNodePositions(settings2, "testMesh", referenceNodePositions);
}

TEST(MeshTest, ReadNodePositionsAreCorrectInlineMesh)
{
  // inline mesh with node positions
  std::string pythonConfig = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [2, 2],
      "nodeDimension": 2,
    "nodePositions": [0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2],
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
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  std::vector<double> referenceNodePositions = {
    0,0,0,  1,0,0,  2,0,0,
    0,1,0,  1,1,0,  2,1,0,
    0,2,0,  1,2,0,  2,2,0,
  };
  
  Mesh::NodePositionsTester::compareNodePositions(settings, "anonymous0", referenceNodePositions);
  
  // inline mesh with automatically generated node positions
  std::string pythonConfig2 = R"(
# Laplace 2D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "physicalExtent": [2.0, 2.0],
      "nElements": [2, 2],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings2(argc, argv, pythonConfig2);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized2(settings2);
  
  Computation computation2(settings2, equationDiscretized2);
  computation2.run();
  
  Mesh::NodePositionsTester::compareNodePositions(settings2, "anonymous0", referenceNodePositions);
  
}

}; // namespace
