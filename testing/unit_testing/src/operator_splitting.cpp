#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "equation/diffusion.h"

// tests only test compilation and do not work inside kdevelop (Fatal error: Could not load dynamic library "lib.so")
/*
TEST(OperatorSplittingTest, Godunov)
{
  std::string pythonConfig = R"(
# Electrophysiology debug
nElements = 5
endTime = 0.1

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58           # membrane capacitance [uF/cm^2]

def setParameters(n_instances, time_step_no, current_time, parameters):
  #print "       > called at ",time_step_no,", time=",current_time, ", n_instances=", n_instances, ", p=",parameters[0]
  
  center_node = n_instances / 2
  
  parameters[0] = 0
  if current_time - int(current_time) < 0.1 and current_time < 10:
    print "parameters len: {}, set I_Stim for node {} to 1200".format(len(parameters),center_node)
    if center_node > 0:
      parameters[center_node-1] = 400.
    parameters[center_node] = 400.
    if center_node < n_instances-1:
      parameters[center_node+1] = 400.
    
config = {
  "Meshes": {
    "MeshFibre": {
      "nElements": nElements,
      "physicalExtend": 5.0,
    },
  },
  "GodunovSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": 1e-1,
    "endTime": endTime,
    "outputData1": True,
    "outputData2": False,

    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out_electro", "outputInterval": 1, "binary": False}
    ],
    "Term1": {      # CellML
      "ExplicitEuler" : {
        "timeStepWidth": 5e-5,
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
        
        "CellML" : {
          "sourceFilename": "cellmlcode.cpp",
          "simdSourceFilename" : "simdcode.cpp",
          #"libraryFilename": "cellml_simd_lib.so",
          
          "numberStates": 57,
          "numberIntermediates": 1,   # intermediates: gamma
          "numberParameters": 2,      # parameters: I_Stim, l_hs
          "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
          "parametersInitialValues": [0.0, 1.0],      # parameters: I_Stim, l_hs
          "meshName": "MeshFibre",
          "prefactor": 1.0,
        },
      },
    },
    "Term2": {     # Diffusion
      "ExplicitEuler" : {
        "timeStepWidth": 1e-5,
        "timeStepOutputInterval": 1e4,
        "FiniteElementMethod" : {
          #"nElements": 0,
          "physicalExtend": 1.0,
          "relativeTolerance": 1e-15,
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
        },
      },
    },
  }
}
)";

  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t), R is from cellml file
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv, pythonConfig);
  
  
  OperatorSplitting::Godunov<
    TimeSteppingScheme::ExplicitEuler<
      CellmlAdapter
    >,
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<1>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<2>,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  >
  problem(settings);
  problem.run();
}


TEST(OperatorSplittingTest, Strang)
{
  std::string pythonConfig = R"(
# Electrophysiology debug
nElements = 5
endTime = 0.1

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58           # membrane capacitance [uF/cm^2]

def setParameters(n_instances, time_step_no, current_time, parameters):
  #print "       > called at ",time_step_no,", time=",current_time, ", n_instances=", n_instances, ", p=",parameters[0]
  
  center_node = n_instances / 2
  
  parameters[0] = 0
  if current_time - int(current_time) < 0.1 and current_time < 10:
    print "parameters len: {}, set I_Stim for node {} to 1200".format(len(parameters),center_node)
    if center_node > 0:
      parameters[center_node-1] = 400.
    parameters[center_node] = 400.
    if center_node < n_instances-1:
      parameters[center_node+1] = 400.
    
config = {
  "Meshes": {
    "MeshFibre": {
      "nElements": nElements,
      "physicalExtend": 5.0,
    },
  },
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-15,
    }
  },
  "StrangSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": 1e-1,
    "endTime": endTime,
    "outputData1": True,
    "outputData2": False,

    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out_electro", "outputInterval": 1, "binary": False}
    ],
    "Term1": {      # CellML
      "ExplicitEuler" : {
        "timeStepWidth": 5e-5,
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
        
        "CellML" : {
          "sourceFilename": "cellmlcode.cpp",
          "simdSourceFilename" : "simdcode.cpp",
          "libraryFilename": "cellml_simd_lib.so",
          
          "numberStates": 57,
          "numberIntermediates": 1,   # intermediates: gamma
          "numberParameters": 2,      # parameters: I_Stim, l_hs
          "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
          "parametersInitialValues": [0.0, 1.0],      # parameters: I_Stim, l_hs
          "meshName": "MeshFibre",
          "prefactor": 1.0,
        },
      },
    },
    "Term2": {     # Diffusion
      "ExplicitEuler" : {
        "timeStepWidth": 1e-5,
        "timeStepOutputInterval": 1e4,
        "FiniteElementMethod" : {
          #"nElements": 0,
          "physicalExtend": 1.0,
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "linearSolver",
        },
      },
    },
  }
}
)";

  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t), R is from cellml file
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv, pythonConfig);
  
  
  OperatorSplitting::Strang<
    TimeSteppingScheme::ExplicitEuler<
      CellmlAdapter
    >,
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<1>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<2>,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  >
  problem(settings);
  problem.run();
}
*/