import sys, os
import timeit
import argparse
import importlib
import distutils.util

# parse rank arguments
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
from helper import *

# meshes

# add neuron meshes
spindle_meshes = {
  "muscleSpindleMesh": {
    # each muscle has the same number of muscle spindels
    "nElements" :         (2*variables.n_muscle_spindles)-1 if n_ranks == 1 else (2*variables.n_muscle_spindles)*n_ranks,
    "physicalExtent":     1, # has no special meaning. only to seperate data points in paraview
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  }
}
variables.meshes.update(spindle_meshes)

# define the config dict
config = {
  "scenarioName":                   variables.scenario_name,    # scenario name which will appear in the log file
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/" + variables.scenario_name + "/mappings_between_meshes.txt",
  "Meshes":                variables.meshes,
  # meshes know their coordinates and the mapping happens automatically. We could also define other parameters such as a mapping tolerance
  # "MappingsBetweenMeshes": {"muscle1_fiber{}".format(f) : ["muscle1Mesh", "muscle1Mesh_quadratic"] for f in range(variables.n_fibers_total)},
  "Solvers": {},


  # connections of the slots, identified by slot name
  # "connectedSlots": [
  #   # global slots only support named slots (connectedSlotsTerm1To2 also allows indices)

  #   # use global slot, because automatic connection of "Razumova/activestress" does not work for some reason
  #   # "Razumova/activestress" from CellML to Muscle contaction solver
  #   ("m1gout", "m1g_in"),
  #   ("m2gout", "m2g_in"),

  #   # lambda and derived values (by MapDofs) -> input of muscle splindel simulation
  #   ("ms0",    "ms_in0"),
  #   ("ms1",    "ms_in1"),
  #   ("ms2",    "ms_in2"),
  #   ("ms3",    "ms_in3"),
  #   ("ms4",    "ms_in4"),

  #   ("in_g",   "in_in"),
  # ],


    # muscle spindles

  #   # mapping muscle spindles output -> motor neuron input
  # "MapDofs": {
  #   "description":                "muscle_spindles_to_motoneurons",   # description that will be shown in solver structure visualization
  #   "nAdditionalFieldVariables":  1,                                  # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
  #   "additionalSlotNames":        ["ms&in"],
  #   "meshName":                   "muscleSpindleMesh",                   # the mesh on which the additional field variables will be defined
  #   "beforeComputation": None,
  #   "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
  #     # each motoneuron gets input from all muscle spindels in the muscle
  #     # ms_out = primary_afferent is max(0, value) an can therefore be zero
  #     {
  #       "fromConnectorSlot":                "ms_out",
  #       "toConnectorSlots":                 ["ms&in"],
  #       "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
  #       "toSlotConnectorArrayIndex":        0,
  #       "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
  #       "fromDofNosNumbering":              "local",
  #       "toDofNosNumbering":                "local",
  #       "dofsMapping":                      None,
  #       "inputDofs":                        muscle1_spindle_indices,   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleMesh"
  #       "outputDofs":                       None, #[in_ms_mesh_muscle1_motoneuron_indices],   # [0,1,...,n_motor_neurons], this is for mesh "motoneuronMesh"
  #       "callback":                         variables.callback_muscle_spindles_to_motoneurons,
  #     }
  #   ],

  "Heun": {
    "description":                  "muscle spindle",
    "timeStepWidth":                variables.dt_muscle_spindles,
    "logTimeStepWidthAsKey":        "dt_muscle_spindles",
    "durationLogKey":               "duration_muscle_spindles",
    "initialValues":                [],
    "timeStepOutputInterval":       500,
    "inputMeshIsGlobal":            True,
    "dirichletBoundaryConditions":  {},
    "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
    "nAdditionalFieldVariables":    0,
    "additionalSlotNames":          [],

    # cellml model of muscle spindle
    "CellML" : {
      "modelFilename":                          variables.muscle_spindle_cellml_file,           # input C++ source file or cellml XML file
      "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
      "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation

      # optimization parameters
      "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
      "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
      "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
      "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.

      # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
      "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
      #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
      "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
      "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
      "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
      "setSpecificStatesRepeatAfterFirstCall":  0.01,                # TODO:                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
      "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
      "additionalArgument":                     None,

      "handleResultFunction":                   None, #handle_result,              # callback function that gets all current values and can do something with them
      "handleResultCallInterval":               1,                          # interval in which handle_result will be called
      "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result

      "mappings":                               variables.muscle_spindle_mappings,              # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
      "parametersInitialValues":                variables.muscle_spindle_parameters_initial_values,  # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]

      "meshName":                               "muscleSpindleMesh",
      "stimulationLogFilename":                 "out/stimulation.log",

      # output writer for states, algebraics and parameters
      "OutputWriter" : [
        {"format": "Paraview",   "outputInterval": int(2./variables.dt_muscle_spindles*variables.output_timestep_spindles), "filename": "out/" + variables.scenario_name + "/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
        {"format": "PythonFile", "outputInterval": int(2./variables.dt_muscle_spindles*variables.output_timestep_spindles), "filename": "out/" + variables.scenario_name + "/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
      ]
    }
  }
  
}


        
        
    
    
  
  

