#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "equation/diffusion.h"
#include "../utility.h"

TEST(DiffusionTest, ExplicitEuler1D)
{
  std::string pythonConfig = R"(

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)


# Diffusion 1D
n = 5
config = {
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
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

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_0000004.py", referenceOutput);

}

TEST(DiffusionTest, Heun1D)
{
  std::string pythonConfig = R"(

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)


# Diffusion 1D
n = 5
config = {
  "Heun" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d_heun", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::Heun<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_heun_0000004.py", referenceOutput);

}

TEST(DiffusionTest, ImplicitEuler1D)
{
  std::string pythonConfig = R"(

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)


# Diffusion 1D
n = 5
config = {
  "ImplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d_implicit", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::ImplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [2.0429559072490386, 2.2518627527228317, 3.8477024200726957, 4.495048744910582, 2.3533552300645306, 2.0611026396733574]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_implicit_0000004.py", referenceOutput);

}

TEST(DiffusionTest, ImplicitEuler1DPOD)
{
  
  std::ofstream file("snapshots.csv");
  file << "2.0,2.0,4.0,5.0,2.0,2.0" << std::endl
    << "2.0034121211861375,2.058007835934967,3.968881360409179,4.883936362778263,2.0849688767570007,2.004998592855786" << std::endl
    << "2.0097405880834756,2.1119059414403125,3.9380119204269333,4.776057992869434,2.1620246733762216,2.0142679983785285" << std::endl
    << "2.018648046962458,2.161994819687798,3.9074790263886396,4.6756844098820585,2.2319744389088783,2.0271049027743073" << std::endl
    << "2.0297760090678523,2.208576359730558,3.8773622943996906,4.582198310705352,2.2955263252435896,2.0429232071836125" << std::endl
    << "2.0428033920520487,2.2519271827902805,3.8477283146699097,4.495039321347134,2.353315231841094,2.061208635718552" << std::endl;
  file.close();
    
  
  std::string pythonConfigReduction = R"(
    # Diffusion 1D POD
n = 5   # number of elements
k = 5   # number of the reduced modes is equal to k+1 because there are 6 snapshots (initial value and five time steps) avalable from the full order model

config = {
  "ModelOrderReduction": {
    "nRowsSnapshots" : n+1,
    "nReducedBases" : k,
      "snapshots" :"./snapshots.csv",
    "ImplicitEuler" : {
       "numberTimeSteps": 5,
       "endTime": 0.1,
       "initialValues": [2,2,4,5,2,2],
       "FiniteElementMethod" : {
          "nElements": n,
          "physicalExtent": 4.0,
          "relativeTolerance": 1e-15,
          "diffusionTensor": [5.0],
       },
       "OutputWriter" : [
         #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
         {"format": "PythonFile", "filename": "diffusion1d_pod_full", "outputInterval": 1, "binary":False}
       ]
    },
    "ImplicitEulerReduced" : {
      "numberTimeSteps": 5,
      "endTime": 0.1,
      "initialValues": [2,2,4,5,2,2],
      "FiniteElementMethod" : {
        "nElements": k,
        "physicalExtent": 4.0,
        "relativeTolerance": 1e-15,
        "diffusionTensor": [5.0],
      },
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
        {"format": "PythonFile", "filename": "diffusion1d_pod_reduced", "outputInterval": 1, "binary":False}
      ]
    },
  },
}
)";

  // problem using reduction
  DihuContext settings(argc, argv, pythonConfigReduction);

  ModelOrderReduction::ImplicitEulerReduced<
    TimeSteppingScheme::ImplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<1>,
        BasisFunction::LagrangeOfOrder<>,
        Quadrature::None,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  > problem(settings);

  problem.run();

  // problem with no reduction
  std::string pythonConfigNoReduction = R"(
# Diffusion 1D
n = 5   # number of elements

config = {
  "ImplicitEuler" : {
     "numberTimeSteps": 5,
     "endTime": 0.1,
     "initialValues": [2,2,4,5,2,2],
     "FiniteElementMethod" : {
        "nElements": n,
        "physicalExtent": 4.0,
        "relativeTolerance": 1e-15,
        "diffusionTensor": [5.0],
     },
     "OutputWriter" : [
       #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
       {"format": "PythonFile", "filename": "diffusion1d_implicit", "outputInterval": 1, "binary":False}
     ]
  }
}
)";
  DihuContext settings2(argc, argv, pythonConfigNoReduction);

  TimeSteppingScheme::ImplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem2(settings2);

  problem2.run();

  // load file contents of reference problem
  std::ifstream outputFile("diffusion1d_implicit_0000004.py");
  std::string referenceOutput((std::istreambuf_iterator<char>(outputFile)), (std::istreambuf_iterator<char>()));

  // compare to full output of POD problem
  assertFileMatchesContent("diffusion1d_pod_full_0000004.py", referenceOutput);
}

TEST(DiffusionTest, CrankNicolson1D)
{
  std::string pythonConfig = R"(
    
import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)
    
    
# Diffusion 1D
n = 5
config = {
  "CrankNicolson" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "relativeTolerance": 1e-15,
    "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d_CN", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::CrankNicolson<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_CN_0000004.py", referenceOutput);

}

TEST(DiffusionTest, ExplicitEuler1DStructuredDeformable)
{
  std::string pythonConfig = R"(

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)


# Diffusion 1D
n = 5
config = {
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<1>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_0000004.py", referenceOutput);

}

TEST(DiffusionTest, Heun1DStructuredDeformable)
{
  std::string pythonConfig = R"(

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)


# Diffusion 1D
n = 5
config = {
  "Heun" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d_heun", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::Heun<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::Gauss<2>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_heun_0000004.py", referenceOutput);

}

TEST(DiffusionTest, ImplicitEuler1DStructuredDeformable)
{
  std::string pythonConfig = R"(

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)


# Diffusion 1D
n = 5
config = {
  "ImplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d_implicit", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::ImplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::Gauss<2>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_implicit_0000004.py", referenceOutput);

}

TEST(DiffusionTest, CrankNicolson1DStructuredDeformable)
{
  std::string pythonConfig = R"(
    
import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)
    
    
# Diffusion 1D
n = 5
config = {
  "CrankNicolson" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "relativeTolerance": 1e-15,
    "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d_implicit", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::CrankNicolson<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::Gauss<2>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredDeformable\", \"dimension\": 1, \"nElementsGlobal\": [5], \"nElementsLocal\": [5], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0, 0.8, 1.6, 2.4000000000000004, 3.2, 4.0]}, {\"name\": \"y\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {\"name\": \"z\", \"values\": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"0\", \"values\": [1.9161287833235827, 2.4114632239553027, 3.842059137806608, 4.19088848006689, 2.6521927547112885, 1.8906640235962393]}]}], \"timeStepNo\": 5, \"currentTime\": 0.1}";
  assertFileMatchesContent("out_diffusion1d_implicit_0000004.py", referenceOutput);

}

TEST(DiffusionTest, Compiles2D)
{
  std::string pythonConfig = R"(
# Diffusion 2D
n = 5

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(n+1) + x
    iv[i] = 1.0

config = {
  "ExplicitEuler" : {
    "initialValues": iv,
    "numberTimeSteps": 5000,
    "endTime": 20.0,
    "FiniteElementMethod" : {
      "nElements": [n,n],
      "physicalExtent": [4.0,4.0],
      "relativeTolerance": 1e-15,
    },
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
