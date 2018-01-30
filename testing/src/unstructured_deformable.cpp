#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

namespace SpatialDiscretization
{
  
TEST(UnstructuredDeformableTest, ReadExfile)
{
#ifdef NDEBUG
  std::cout<< "NDEBUG is defined"<<std::endl;
#else
  std::cout<< "NDEBUG is not defined"<<std::endl;
#endif
  
  std::string pythonConfig = R"(
# Laplace 3D
    
# boundary conditions
bc = {}
bc[0] = 1.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "exelem": "left_biceps_brachii.exelem",
    "exnode": "left_biceps_brachii.exnode",
    "OutputWriter" : [
      {"format": "Exfile", "interval": 1, "filename": "out"},
      #{"format": "Python", "filename": "p"}
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::UnstructuredDeformable<3>,
    BasisFunction::Hermite,
    Integrator::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
}

TEST(UnstructuredDeformableTest, CreateFromConfig2D)
{
  std::string pythonConfig = R"(
# Laplace 2D
    
# boundary conditions
bc = {}
bc[0] = 1.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "nodePositions": [[0,0,0], [1,0], [2,0,0], [0,1], [1,1], [2,1], [0,2], [1,2], [2,2]],  # 3x3 nodes, 4 elements
    "elements": [[[0,0], 1, [3,0], 4], [[1,0], [2,0], [4,0], [5,0]], [[3,0], [4,0], [6,0], [7,0]], [[4,0], [5,0], [7,0], [8,0]]],   # each node is [node no, version-at-that-node no] or just node-no then it assumes version no 0
    "OutputWriter" : [
      {"format": "Exfile", "interval": 1, "filename": "init_from_config"},
    ]
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::UnstructuredDeformable<2>,
    BasisFunction::Hermite,    // 4 nodes per element
    Integrator::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
}




};

