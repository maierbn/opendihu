#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "equation/diffusion.h"

TEST(DiffusionTest, Compiles1D)
{
  std::string pythonConfig = R"(
# Diffusion 1D
n = 5
config = {
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 1.0,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);
  
  problem.run();
}

TEST(DiffusionTest, Compiles2D)
{
  std::string pythonConfig = R"(
# Diffusion 2D
n = 50

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(n+1) + x
    iv[i] = 1.0

config = {
  "FiniteElementMethod" : {
    "nElements": [n,n],
    "physicalExtend": [4.0,4.0],
    "relativeTolerance": 1e-15,
  },
  "ExplicitEuler" : {
    "initialValues": iv,
    "numberTimeSteps": 5000,
    "endTime": 20.0,
  },
  "OutputWriter" : [
    {"format": "PythonFile", "filename": "out_diffusion2d", "frequency": 100, "binary": False}
  ]
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);
}

