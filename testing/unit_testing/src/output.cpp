#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "utility.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

namespace SpatialDiscretization
{
  
TEST(OutputTest, UnstructuredDeformable)
{
  std::string pythonConfig = R"(
# Laplace 2D
n=4

def callback(config):
  with open("result_callback","w") as f:
    f.write(str(config))
    
  import numpy as np
  a = np.load("out_binary.npy")
  with open("result_binary","w") as f:
    f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "DirichletBoundaryCondition": {0:1.0},
    "relativeTolerance": 1e-15,
    "nodePositions": [[0,0,0], [1,0], [2,0,0], [0,1], [1,1], [2,1], [0,2], [1,2], [2,2]],  # 3x3 nodes, 4 elements
    "elements": [[0, 1, 3, 4], [1, 2, 4, 5], [3, 4, 6, 7], [4, 5, 7, 8]],   # each node is [node no, version-at-that-node no] or just node-no then it assumes version no 0
    "OutputWriter" : [
      {"format": "PythonFile", "filename" : "out_binary", "binary" : True},
      {"format": "PythonFile", "filename" : "out_txt", "binary" : False},
      {"format": "PythonCallback", "callback": callback},
    ]
  }
}
)";
  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  std::string referenceOutput = "{'timeStepNo': -1, 'basisFunction': 'Lagrange', 'meshType': 'UnstructuredDeformable', 'nElements': 4, 'currentTime': 0.0, 'basisOrder': 1, 'onlyNodalValues': True, 'data': [{'name': 'geometry', 'components': [{'values': [0.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 2.0], 'name': 'x'}, {'values': [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0], 'name': 'y'}, {'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'name': 'z'}]}, {'name': 'solution', 'components': [{'values': [0.9999999999999998, 0.9999999999999991, 0.9999999999999991, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999989], 'name': '0'}]}, {'name': 'rhs', 'components': [{'values': [1.0, -0.16666666666666669, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0, 0.0], 'name': '0'}]}], 'dimension': 2}";
  std::string referenceOutputSolution = "[ 1.  1.  1.  1.  1.  1.  1.  1.  1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput);
  assertFileMatchesContent("out_txt.py", referenceOutput);
  assertFileMatchesContent("result_binary", referenceOutputSolution);
}
 
TEST(OutputTest, StructuredDeformable)
{
  std::string pythonConfig = R"(
# Laplace 2D
n=4

def callback(config):
  print "config=",config
  with open("result_callback","w") as f:
    f.write(str(config))
    
  import numpy as np
  a = np.load("out_binary.npy")
  print "binary=",a
  with open("result_binary","w") as f:
    f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "DirichletBoundaryCondition": {0:1.0},
    "relativeTolerance": 1e-15,
    "nElements": [2,2],
    "nodeDimension": 2,
    "nodePositions": [0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2],  # 3x3 nodes, 4 elements
    "OutputWriter" : [
      {"format": "PythonFile", "filename" : "out_binary", "binary" : True},
      {"format": "PythonFile", "filename" : "out_txt", "binary" : False},
      {"format": "PythonCallback", "callback": callback},
    ]
  }
}
)";
  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  std::string referenceOutput = "{'timeStepNo': -1, 'basisFunction': 'Lagrange', 'meshType': 'StructuredDeformable', 'nElements': [2, 2], 'currentTime': 0.0, 'basisOrder': 1, 'onlyNodalValues': True, 'data': [{'name': 'geometry', 'components': [{'values': [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0], 'name': 'x'}, {'values': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0], 'name': 'y'}, {'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'name': 'z'}]}, {'name': 'solution', 'components': [{'values': [1.0, 0.999999999999999, 0.9999999999999987, 0.9999999999999991, 0.9999999999999988, 0.9999999999999987, 0.9999999999999984, 0.9999999999999984, 0.9999999999999989], 'name': '0'}]}, {'name': 'rhs', 'components': [{'values': [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0], 'name': '0'}]}], 'dimension': 2}";
  std::string referenceOutputSolution = "[ 1.  1.  1.  1.  1.  1.  1.  1.  1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput);
  assertFileMatchesContent("out_txt.py", referenceOutput);
  assertFileMatchesContent("result_binary", referenceOutputSolution);
}

// segfault, python memory problem in RegularFixed
/*
TEST(OutputTest, RegularFixed)
{
  std::string pythonConfig = R"(
# Laplace 2D
n=4

def callback(config):
  print "config=",config
  with open("result_callback","w") as f:
    f.write(str(config))
    
  import numpy as np
  a = np.load("out_binary.npy")
  print "binary=",a
  with open("result_binary","w") as f:
    f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": [4, 4],
    "physicalExtent": [4.0, 4.0],
    "DirichletBoundaryCondition": {0:1.0},
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "PythonFile", "filename" : "out_binary", "binary" : True},
      {"format": "PythonFile", "filename" : "out_txt", "binary" : False},
      {"format": "PythonCallback", "callback": callback},
    ]
  }
}
)";
  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Static::Laplace
  > equationDiscretized(settings);
  
  Computation computation(settings, equationDiscretized);
  computation.run();
  
  std::string referenceOutput = "{'timeStepNo': -1, 'currentTime': 0.0, 'nElements': [4, 4], 'meshType': 'RegularFixed', 'data': [{'name': 'geometry', 'components': [{'values': [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0], 'name': 'x'}, {'values': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0], 'name': 'y'}, {'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'name': 'z'}]}, {'name': 'solution', 'components': [{'values': [1.0, 0.9999999999999989, 0.9999999999999987, 0.9999999999999989, 0.9999999999999989, 0.9999999999999987, 0.9999999999999987, 0.9999999999999989, 0.9999999999999987], 'name': '0'}]}, {'name': 'rhs', 'components': [{'values': [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0], 'name': '0'}]}], 'dimension': 2}";
  std::string referenceOutputSolution = "[ 1.  1.  1.  1.  1.  1.  1.  1.  1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput);
  assertFileMatchesContent("out_txt.py", referenceOutput);
  assertFileMatchesContent("result_binary", referenceOutputSolution);
}
*/
};

