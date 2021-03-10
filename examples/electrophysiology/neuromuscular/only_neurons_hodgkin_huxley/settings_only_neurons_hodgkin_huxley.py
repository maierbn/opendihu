# sensory organs, interneurons, motoneurons, multidomain mockup (no electrophyisiology), mechanics
#
# [muscle spindle model] -> [motor neuron model] <- mapping by callback <- [interneuron model] <- mapping by callback <- [golgi tendon organ model]
#                            +-> electro-mechanics model

import numpy as np
import scipy.stats
import pickle
import sys,os 
import timeit
import argparse
import importlib
import struct
sys.path.insert(0, "..")

# own MPI rank no and number of MPI ranks
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there  
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

# if first argument contains "*.py", it is a custom variable definition file, load these values
if ".py" in sys.argv[0]:
  variables_path_and_filename = sys.argv[0]
  variables_path,variables_filename = os.path.split(variables_path_and_filename)  # get path and filename 
  sys.path.insert(0, os.path.join(script_path,variables_path))                    # add the directory of the variables file to python path
  variables_module,_ = os.path.splitext(variables_filename)                       # remove the ".py" extension to get the name of the module
  
  if rank_no == 0:
    print("Loading variables from \"{}\".".format(variables_path_and_filename))
    
  custom_variables = importlib.import_module(variables_module, package=variables_filename)    # import variables module
  variables.__dict__.update(custom_variables.__dict__)
  sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed
else:
  if rank_no == 0:
    print("Warning: There is no variables file, e.g:\n ./only_neurons ../settings_only_neurons_hodgkin_huxley.py normal.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='only_neurons')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

# output information of run
if rank_no == 0:
  print("scenario_name: {},   n_ranks: {},  end_time: {}".format(variables.scenario_name, n_ranks, variables.end_time))
  print("dt_neuron_transfer:  {:0.0e}".format(variables.dt_neuron_transfer))
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# add neuron meshes
variables.meshes = {
  "motoneuronMesh": {
    "nElements" :         variables.n_motoneurons-1 if n_ranks == 1 else variables.n_motoneurons*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "motoneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "interneuronMesh": {
    "nElements" :         variables.n_interneurons-1 if n_ranks == 1 else variables.n_interneurons*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "interneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleMesh": {
    "nElements" :         variables.n_muscle_spindles-1 if n_ranks == 1 else variables.n_muscle_spindles*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleAndInterneuronMesh": {
    "nElements" :         variables.n_muscle_spindles+variables.n_interneurons-1 if n_ranks == 1 else (variables.n_muscle_spindles+variables.n_interneurons)*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle_interneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "golgiTendonOrganMesh": {
    "nElements" :         variables.n_golgi_tendon_organs-1 if n_ranks == 1 else variables.n_golgi_tendon_organs*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "golgi_tendon_organ",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
}
  
config = {
  "scenarioName":                   variables.scenario_name,              # scenario name which will be included in the log file logs/log.csv, to identify the current run
  "logFormat":                      "csv",                                # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",               # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes_log.txt",    # log file for mappings 
  "Meshes":                         variables.meshes,                     # all used meshes
  
  # connections of the slots, identified by slot name
  "connectedSlots": [
    ("mn_out", "mn"),
    ("in_g",   "in_in"),
    ("msin_s", "msin_i"),
    ("msin_i", "msin_m"),
    ("gt",     "gt_in"),
    ("ms",     "ms_in"),
  ],
  
  # global coupling scheme
  "MultipleCoupling": {
    "endTime":                variables.end_time,                         # end time of the simulation
    "timeStepWidth":          variables.dt_neuron_transfer,               # time step width of the data transfer between the sub solvers
    "logTimeStepWidthAsKey":  "dt_neuron_transfer",                       # string under which the timestep width will be stored in the log file
    "durationLogKey":         "duration_total",                           # string under which the total duration will be stored in the log file
    "timeStepOutputInterval": 500,                                        # how often to print the current time step
    "connectedSlotsTerm1To2": None,                                       # this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead
    "connectedSlotsTerm2To1": None,                                       # this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead

    # muscle spindles, golgi tendon organs and interneurons
    "Term1": {
      # mapping muscle spindles output -> motor neuron signals
      "MapDofs": {
        "description":                "muscle_spindles_to_motoneurons",   # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  1,                                  # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["msin_s"],
        "meshName":                   "muscleSpindleAndInterneuronMesh",  # the mesh on which the additional field variables will be defined
        "beforeComputation": None,
        "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                "ms_out",
            "toConnectorSlots":                 ["msin_s"],
            "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
            "toSlotConnectorArrayIndex":        0,
            "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "local",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        list(range(variables.n_muscle_spindles)),   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleMesh"
            "outputDofs":                       [list(range(variables.n_muscle_spindles))],   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleAndInterneuronMesh"
            "callback":                         variables.callback_muscle_spindles_to_motoneurons,
            #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
            #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          }
        ],
               
        # muscle spindle model solver
        "Heun" : {
          "description":                  "muscle spindle",
          "timeStepWidth":                variables.dt_muscle_spindles,
          "logTimeStepWidthAsKey":        "dt_muscle_spindles",
          "durationLogKey":               "duration_muscle_spindles",
          "initialValues":                [],
          "timeStepOutputInterval":       1e4,
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
            "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
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
              {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_muscle_spindles*variables.output_timestep_neurons), "filename": "out/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
            ]
          }
        }
      }
    },
    
    "Term2": {
      # mapping Golgi tendon organs -> interneurons
      "MapDofs": {
        "description":                "golgi_tendon_o._to_interneurons",       # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  1,                             # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["in_g"],
        "meshName":                   "interneuronMesh",             # the mesh on which the additional field variables will be defined
        "beforeComputation": None,
        "afterComputation": [                                        # transfer/mapping of dofs that will be performed after the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                "gt_out",
            "toConnectorSlots":                 "in_g",
            "fromSlotConnectorArrayIndex":      0,                   # which fiber/compartment
            "toSlotConnectorArrayIndex":        0,
            "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "local",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        list(range(variables.n_golgi_tendon_organs)),   # [0,1,...,n_golgi_tendon_organs], this is for mesh "golgiTendonOrganMesh"
            "outputDofs":                       [list(range(variables.n_interneurons))],          # [0,1,...,n_interneurons]         this is for mesh "interneuronMesh"
            "callback":                         variables.callback_golgi_tendon_organs_to_interneurons,
            #"thresholdValue":                   20,                 # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
            #"valueToSet":                       20,                 # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          }
        ],
        # Golgi tendon organ model solver
        "Heun" : {
          "description":                  "Golgi tendon organs",
          "timeStepWidth":                variables.dt_golgi_tendon_organs,
          "logTimeStepWidthAsKey":        "dt_golgi_tendon_organs",
          "durationLogKey":               "duration_golgi_tendon_organ",
          "initialValues":                [],
          "timeStepOutputInterval":       1e4,
          "inputMeshIsGlobal":            True,
          "dirichletBoundaryConditions":  {},
          "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
          "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
          "nAdditionalFieldVariables":    0,
          "additionalSlotNames":          [],
              
          # cellml model of golgi tendon organs
          "CellML" : {
            "modelFilename":                          variables.golgi_tendon_organ_cellml_file,       # input C++ source file or cellml XML file
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
            "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
            "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
            "additionalArgument":                     None,
            
            "mappings":                               variables.golgi_tendon_organ_mappings,          # mappings between parameters and algebraics/constants and between connectorSlots and states, algebraics or parameters, they are defined in helper.py
            "parametersInitialValues":                variables.golgi_tendon_organ_parameters_initial_values,    # # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
            
            "meshName":                               "golgiTendonOrganMesh",                         
            "stimulationLogFilename":                 "out/stimulation.log",

            # output writer for states, algebraics and parameters                
            "OutputWriter" : [
              {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_golgi_tendon_organs*variables.output_timestep_neurons), "filename": "out/golgi_tendon_organs", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
            ]
          }
        }
      }
    },
    
    "Term3": {
      # mapping interneurons -> input for motor neurons
      "MapDofs": {
        "description":                "interneurons_to_motoneurons",  # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["msin_i"],
        "meshName":                   "muscleSpindleAndInterneuronMesh", # the mesh on which the additional field variables will be defined
        "beforeComputation": None,                                    # transfer/mapping of dofs that will be performed before the computation of the nested solver
        "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                "in_out",
            "toConnectorSlots":                 "msin_i",
            "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
            "toSlotConnectorArrayIndex":        0,
            "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "local",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        list(range(variables.n_interneurons)),   # [0,1,...,n_interneurons]
            "outputDofs":                       [list(range(variables.n_muscle_spindles,variables.n_muscle_spindles+variables.n_interneurons))],          # [n_muscle_spindles,n_muscle_spindles+1,...,n_muscle_spindles+n_interneurons]
            "callback":                         variables.callback_interneurons_to_motoneurons,
            #"thresholdValue":                   20,                    # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
            #"valueToSet":                       20,                   # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          }
        ],
               
        # interneuron solver
        "Heun" : {
          "description":                  "interneurons",
          "timeStepWidth":                variables.dt_interneuron,
          "logTimeStepWidthAsKey":        "dt_interneuron",
          "durationLogKey":               "duration_interneuron",
          "initialValues":                [],
          "timeStepOutputInterval":       1e4,
          "inputMeshIsGlobal":            True,
          "dirichletBoundaryConditions":  {},
          "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
          "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
          "nAdditionalFieldVariables":    0,
          "additionalSlotNames":          [],
              
          # cellml model of interneurons
          "CellML" : {
            "modelFilename":                          variables.interneuron_cellml_file,              # input C++ source file or cellml XML file
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
            "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
            "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
            "additionalArgument":                     None,
            
            "mappings":                               variables.interneuron_mappings,                 # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
            "parametersInitialValues":                variables.interneuron_parameters_initial_values, # # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
            
            "meshName":                               "interneuronMesh",                           
            "stimulationLogFilename":                 "out/stimulation.log",

            # output writer for states, algebraics and parameters                
            "OutputWriter" : [
              {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_interneuron*variables.output_timestep_neurons), "filename": "out/interneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
            ]
          }
        }
      }
    },
    
    "Term4": {
      # mapping motor neuron signals + cortical input to actual inputs
      "MapDofs": {
        "description":                "motoneurons_input",   # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["msin_m"],
        "meshName":                   "muscleSpindleAndInterneuronMesh",               # the mesh on which the additional field variables will be defined
        "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                "msin_m",
            "toConnectorSlots":                 "mn_in",
            "fromSlotConnectorArrayIndex":      0,
            "toSlotConnectorArrayIndex":        0,
            "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "local",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        list(range(variables.n_muscle_spindles+variables.n_interneurons)),   # [0,1,...,n_muscle_spindles,...,n_muscle_spindles+n_interneurons]
            "outputDofs":                       [list(range(variables.n_motoneurons))],   # [0,1,...,n_motoneurons]
            "callback":                         variables.callback_motoneurons_input,
            #"thresholdValue":                   20,                   # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
            #"valueToSet":                       20,                   # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          }
        ],
        "afterComputation": None,
                        
        # motoneuron solver
        "Heun" : {
          "description":                  "motoneurons",
          "timeStepWidth":                variables.dt_motoneuron,
          "logTimeStepWidthAsKey":        "dt_motoneuron",
          "durationLogKey":               "duration_motoneuron",
          "initialValues":                [],
          "timeStepOutputInterval":       1e4,
          "inputMeshIsGlobal":            True,
          "dirichletBoundaryConditions":  {},
          "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
          "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
          "nAdditionalFieldVariables":    0,
          "additionalSlotNames":          [],
              
          # cellml model of motorneuron
          "CellML" : {
            "modelFilename":                          variables.motoneuron_cellml_file,               # input C++ source file or cellml XML file
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
            "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
            "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
            "additionalArgument":                     None,
            
            "mappings":                               variables.motoneuron_mappings,                  # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
            "parametersInitialValues":                variables.motoneuron_parameters_initial_values, # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
            
            "meshName":                               "motoneuronMesh",                               # use the linear mesh, it was partitioned by the helper.py script which called opendihu/scripts/create_partitioned_meshes_for_settings.py
            "stimulationLogFilename":                 "out/stimulation.log",

            # output writer for states, algebraics and parameters                
            "OutputWriter" : [
              {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "filename": "out/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
            ]
          }
        }
      }
    },
        
    # electrophysiology and mechanics solver
    "Term5": {
      # map from λ in the 3D mesh to muscle spindles input
      "MapDofs": {
        "description":                "muscle_spindles_input", # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["ms"],
        "meshName":                   "muscleSpindleMesh",            # the mesh on which the additional field variables will be defined
        "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                0,
            "toConnectorSlots":                 "ms",
            "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment, this does not matter here because all compartment meshes have the same displacements
            "toSlotConnectorArrayIndex":        0,
            "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "global",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        None,          # nodes are set at bottom in helper.py
            "outputDofs":                       [list(range(variables.n_muscle_spindles))],   # [0,1,...,n_muscle_spindles]
            "callback":                         variables.callback_muscle_spindles_input,
            #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
            #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          }
        ],
        "afterComputation":  None,
          
        # map from λ in the 3D mesh to golgi tendon organs
        "MapDofs": {
          "description":                "golgi_tendon_organs_input",      # description that will be shown in solver structure visualization
          "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
          "additionalSlotNames":        ["gt"],
          "meshName":                   "golgiTendonOrganMesh",               # the mesh on which the additional field variables will be defined
          "beforeComputation":          None, 
          "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
            {                                                 
              "fromConnectorSlot":                0,
              "toConnectorSlots":                 "gt",
              "fromSlotConnectorArrayIndex":      0,                   # which fiber/compartment, this does not matter here because all compartment meshes have the same displacements
              "toSlotConnectorArrayIndex":        0,
              "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
              "fromDofNosNumbering":              "global",
              "toDofNosNumbering":                "local",
              "dofsMapping":                      None,
              "inputDofs":                        None,          # nodes are set at bottom in helper.py
              "outputDofs":                       [list(range(variables.n_golgi_tendon_organs))],   # [0,1,...,n_golgi_tendon_organs]
              "callback":                         variables.callback_golgi_tendon_organs_input,
              #"thresholdValue":                   20,                 # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
              #"valueToSet":                       20,                 # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
            }
          ],
              
          # map from motoneuronMesh to stimulated nodes, this is not used here because there is no muscle mesh
          "MapDofs": {
            "description":                "motoneurons->stimulated nodes",  # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["mn"],
            "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
            "beforeComputation":          None, 
            "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
              {
                "fromConnectorSlot":                0, #2,
                "toConnectorSlots":                 0,
                "fromSlotConnectorArrayIndex":      0,
                "toSlotConnectorArrayIndex":        compartment_no,      # which motor unit
                "mode":                             "localSetIfAboveThreshold",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "global",
                "dofsMapping":                      None,
                "inputDofs":                        compartment_no,
                "outputDofs":                       None,
                #"callback":                         variables.callback_motoneuron_output,
                "thresholdValue":                   20,                    # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                "valueToSet":                       20,                   # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
            for compartment_no in range(variables.n_compartments)],
              
            # electro-mechanics solver
            # Dummy
          }
        }
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))

