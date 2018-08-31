#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "../utility.h"
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
    f.write(str(config[0]))
    
  #import numpy as np
  #a = np.load("out_binary_solution.npy")
  #with open("result_binary","w") as f:
  #  f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "initialValues": [0],
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
  
   
  equationDiscretized.run();
  
  std::string referenceOutput = "{'meshType': 'UnstructuredDeformable', 'dimension': 2, 'nElements': 4, 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 2.0]}, {'name': 'y', 'values': [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [0.9999999999999998, 0.9999999999999991, 0.9999999999999991, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999989]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666669, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], 'elementalDofs': [[0, 1, 2, 3], [1, 4, 3, 5], [2, 3, 6, 7], [3, 5, 7, 8]], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput2 = "{\"meshType\": \"UnstructuredDeformable\", \"dimension\": 2, \"nElements\": 4, \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.9999999999999998, 0.9999999999999991, 0.9999999999999991, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999988, 0.9999999999999989]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666669, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"elementalDofs\": [[0, 1, 2, 3], [1, 4, 3, 5], [2, 3, 6, 7], [3, 5, 7, 8]], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput3 = "{'meshType': 'UnstructuredDeformable', 'dimension': 2, 'nElements': 4, 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 2.0]}, {'name': 'y', 'values': [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [0.9999999999999998, 0.9999999999999987, 0.9999999999999988, 0.9999999999999987, 0.9999999999999987, 0.9999999999999989, 0.9999999999999987, 0.9999999999999986, 0.9999999999999986]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666669, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], 'elementalDofs': [[0, 1, 2, 3], [1, 4, 3, 5], [2, 3, 6, 7], [3, 5, 7, 8]], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput4 = "{\"meshType\": \"UnstructuredDeformable\", \"dimension\": 2, \"nElements\": 4, \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [0.9999999999999998, 0.9999999999999987, 0.9999999999999988, 0.9999999999999987, 0.9999999999999987, 0.9999999999999989, 0.9999999999999987, 0.9999999999999986, 0.9999999999999986]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666669, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"elementalDofs\": [[0, 1, 2, 3], [1, 4, 3, 5], [2, 3, 6, 7], [3, 5, 7, 8]], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  //std::string referenceOutputSolution = "[1. 1. 1. 1. 1. 1. 1. 1. 1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput, referenceOutput3);
  assertFileMatchesContent("out_txt.py", referenceOutput2, referenceOutput4);
  //assertFileMatchesContent("result_binary", referenceOutputSolution);
}

TEST(OutputTest, StructuredDeformable)
{
  std::string pythonConfig = R"(
# Laplace 2D

def callback(config):
  print("config=",config)
  with open("result_callback","w") as f:
    f.write(str(config[0]))
    
  #import numpy as np
  #a = np.load("out_binary_solution.npy")
  #print("binary=",a)
  #with open("result_binary","w") as f:
  #  f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "DirichletBoundaryCondition": {0:1.0},
    "initialValues": [0],
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
  
   
  equationDiscretized.run();
  
  std::string referenceOutput = "{'meshType': 'StructuredDeformable', 'dimension': 2, 'nElementsGlobal': [2, 2], 'nElementsLocal': [2, 2], 'beginNodeGlobalNatural': [0, 0], 'hasFullNumberOfNodes': [True, True], 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {'name': 'y', 'values': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [1.0, 0.999999999999999, 0.9999999999999987, 0.9999999999999991, 0.9999999999999988, 0.9999999999999987, 0.9999999999999984, 0.9999999999999984, 0.9999999999999989]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput2 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [2, 2], \"nElementsLocal\": [2, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0, 0.999999999999999, 0.9999999999999987, 0.9999999999999991, 0.9999999999999988, 0.9999999999999987, 0.9999999999999984, 0.9999999999999984, 0.9999999999999989]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput3 = "{'meshType': 'StructuredDeformable', 'dimension': 2, 'nElementsGlobal': [2, 2], 'nElementsLocal': [2, 2], 'beginNodeGlobalNatural': [0, 0], 'hasFullNumberOfNodes': [True, True], 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {'name': 'y', 'values': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [1.0, 0.9999999999999962, 0.9999999999999958, 0.9999999999999964, 0.9999999999999961, 0.9999999999999958, 0.9999999999999957, 0.9999999999999956, 0.9999999999999958]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput4 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [2, 2], \"nElementsLocal\": [2, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0, 0.9999999999999962, 0.9999999999999958, 0.9999999999999964, 0.9999999999999961, 0.9999999999999958, 0.9999999999999957, 0.9999999999999956, 0.9999999999999958]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  //std::string referenceOutputSolution = "[1. 1. 1. 1. 1. 1. 1. 1. 1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput, referenceOutput3);
  assertFileMatchesContent("out_txt.py", referenceOutput2, referenceOutput4);
  //assertFileMatchesContent("result_binary", referenceOutputSolution);
}

TEST(OutputTest, StructuredDeformable2)
{
  std::string pythonConfig = R"(
# Laplace 2D

def callback(config):
  print("config=",config)
  with open("result_callback","w") as f:
    f.write(str(config[0]))
    
  #import numpy as np
  #a = np.load("out_binary_solution.npy")
  #print("binary=",a)
  #with open("result_binary","w") as f:
  #  f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "DirichletBoundaryCondition": {0:1.0},
    "initialValues": [0],
    "relativeTolerance": 1e-15,
    "nElements": [2,2],
    "nodeDimension": 2,
    "nodePositions": [[0,0], 1, [2], [0,1], [1,1,0], [2,1,0], [0,2], [1,2], [2,2]],  # 3x3 nodes, 4 elements
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
  
   
  equationDiscretized.run();
  
  std::string referenceOutput = "{'meshType': 'StructuredDeformable', 'dimension': 2, 'nElementsGlobal': [2, 2], 'nElementsLocal': [2, 2], 'beginNodeGlobalNatural': [0, 0], 'hasFullNumberOfNodes': [True, True], 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {'name': 'y', 'values': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [1.0, 0.999999999999999, 0.9999999999999987, 0.9999999999999991, 0.9999999999999988, 0.9999999999999987, 0.9999999999999984, 0.9999999999999984, 0.9999999999999989]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput2 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [2, 2], \"nElementsLocal\": [2, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0, 0.999999999999999, 0.9999999999999987, 0.9999999999999991, 0.9999999999999988, 0.9999999999999987, 0.9999999999999984, 0.9999999999999984, 0.9999999999999989]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  std::string referenceOutput3 = "{'meshType': 'StructuredDeformable', 'dimension': 2, 'nElementsGlobal': [2, 2], 'nElementsLocal': [2, 2], 'beginNodeGlobalNatural': [0, 0], 'hasFullNumberOfNodes': [True, True], 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {'name': 'y', 'values': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [1.0, 0.9999999999999962, 0.9999999999999958, 0.9999999999999964, 0.9999999999999961, 0.9999999999999958, 0.9999999999999957, 0.9999999999999956, 0.9999999999999958]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput4 = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 2, \"nElementsGlobal\": [2, 2], \"nElementsLocal\": [2, 2], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0, 0.9999999999999962, 0.9999999999999958, 0.9999999999999964, 0.9999999999999961, 0.9999999999999958, 0.9999999999999957, 0.9999999999999956, 0.9999999999999958]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666669, 0.0, -0.16666666666666669, -0.33333333333333337, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  //std::string referenceOutputSolution = "[1. 1. 1. 1. 1. 1. 1. 1. 1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput, referenceOutput3);
  assertFileMatchesContent("out_txt.py", referenceOutput2, referenceOutput4);
  //assertFileMatchesContent("result_binary", referenceOutputSolution);
}


TEST(OutputTest, RegularFixed)
{
  std::string pythonConfig = R"(
# Laplace 2D

def callback(config):
  print("config=",config)
  with open("result_callback","w") as f:
    f.write(str(config[0]))
    
  #import numpy as np
  #a = np.load("out_binary_solution.npy")
  #print("binary=",a)
  #with open("result_binary","w") as f:
  #  f.write(str(a))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": [4, 4],
    "physicalExtent": [4.0, 4.0],
    "initialValues": [0],
    "DirichletBoundaryCondition": {0:1.0},
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "PythonFile", "filename" : "out_binary", "binary" : True},
      {"format": "PythonCallback", "callback": callback},
      {"format": "PythonFile", "filename" : "out_txt", "binary" : False},
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
  
   
  equationDiscretized.run();
  
  std::string referenceOutput = "{'meshType': 'StructuredRegularFixed', 'dimension': 2, 'nElementsGlobal': [4, 4], 'nElementsLocal': [4, 4], 'beginNodeGlobalNatural': [0, 0], 'hasFullNumberOfNodes': [True, True], 'basisFunction': 'Lagrange', 'basisOrder': 1, 'onlyNodalValues': True, 'nRanks': 1, 'ownRankNo': 0, 'data': [{'name': 'geometry', 'components': [{'name': 'x', 'values': [0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0]}, {'name': 'y', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0]}, {'name': 'z', 'values': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {'name': 'solution', 'components': [{'name': '0', 'values': [1.0000000000000004, 1.0000000000000018, 1.0000000000000027, 1.0000000000000036, 1.0000000000000033, 1.0000000000000018, 1.0000000000000022, 1.0000000000000027, 1.0000000000000029, 1.000000000000003, 1.0000000000000027, 1.0000000000000027, 1.0000000000000036, 1.0000000000000033, 1.000000000000003, 1.0000000000000036, 1.000000000000003, 1.0000000000000033, 1.0000000000000027, 1.000000000000003, 1.0000000000000033, 1.0000000000000036, 1.0000000000000033, 1.000000000000003, 1.000000000000003]}]}, {'name': 'rhs', 'components': [{'name': '0', 'values': [1.0, -0.16666666666666666, 0.0, 0.0, 0.0, -0.16666666666666666, -0.3333333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], 'timeStepNo': -1, 'currentTime': 0.0}";
  std::string referenceOutput2 = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 2, \"nElementsGlobal\": [4, 4], \"nElementsLocal\": [4, 4], \"beginNodeGlobalNatural\": [0, 0], \"hasFullNumberOfNodes\": [true, true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.0000000000000004, 1.0000000000000018, 1.0000000000000027, 1.0000000000000036, 1.0000000000000033, 1.0000000000000018, 1.0000000000000022, 1.0000000000000027, 1.0000000000000029, 1.000000000000003, 1.0000000000000027, 1.0000000000000027, 1.0000000000000036, 1.0000000000000033, 1.000000000000003, 1.0000000000000036, 1.000000000000003, 1.0000000000000033, 1.0000000000000027, 1.000000000000003, 1.0000000000000033, 1.0000000000000036, 1.0000000000000033, 1.000000000000003, 1.000000000000003]}]}, {\"name\": \"rhs\", \"components\": [{\"name\": \"0\", \"values\": [1.0, -0.16666666666666666, 0.0, 0.0, 0.0, -0.16666666666666666, -0.3333333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}], \"timeStepNo\": -1, \"currentTime\": 0.0}";
  //std::string referenceOutputSolution = "[1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.\n 1.]";
  
  assertFileMatchesContent("result_callback", referenceOutput);
  assertFileMatchesContent("out_txt.py", referenceOutput2);
  //assertFileMatchesContent("result_binary", referenceOutputSolution);
}

};

