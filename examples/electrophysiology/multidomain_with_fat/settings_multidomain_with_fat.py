# Multicompartment 3D, biceps
#

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
  variables_file = sys.argv[0]
  variables_module = variables_file[0:variables_file.find(".py")]
  
  if rank_no == 0:
    print("Loading variables from {}.".format(variables_file))
    
  custom_variables = importlib.import_module(variables_module)
  variables.__dict__.update(custom_variables.__dict__)
  sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed
else:
  if rank_no == 0:
    print("Warning: There is no variables file, e.g:\n ./multidomain_with_fat ../settings_multidomain_with_fat.py debug.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='multidomain_with_fat')
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)

# parse command line arguments and assign values to variables module
args = parser.parse_args(args=sys.argv[:-2], namespace=variables)

# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]
  
variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:
  
  # create all possible partitionings to the given number of ranks
  optimal_value = n_ranks**(1/3)
  possible_partitionings = []
  for i in range(1,n_ranks):
    for j in [1]:
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = (int)(n_ranks / (i*j))
        performance = (k-optimal_value)**2 + (j-optimal_value)**2 + 1.1*(i-optimal_value)**2
        possible_partitionings.append([i,j,k,performance])
        
  # if no possible partitioning was found
  if len(possible_partitionings) == 0:
    if rank_no == 0:
      print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {} and no automatic partitioning could be done.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    quit()
    
  # select the partitioning with the lowest value of performance which is the best
  lowest_performance = possible_partitionings[0][3]+1
  for i in range(len(possible_partitionings)):
    if possible_partitionings[i][3] < lowest_performance:
      lowest_performance = possible_partitionings[i][3]
      variables.n_subdomains_x = possible_partitionings[i][0]
      variables.n_subdomains_y = possible_partitionings[i][1]
      variables.n_subdomains_z = possible_partitionings[i][2]

# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.0e}".format(variables.dt_0D))
  print("dt_splitting:    {:0.0e}".format(variables.dt_splitting))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

# settings for the multidomain solver
multidomain_solver = {
  "nCompartments":                    variables.n_compartments,           # number of compartments
  "am":                               [variables.get_am(mu_no) for mu_no in range(variables.n_compartments)],   # Am parameter for every motor unit (ration of surface to volume of fibers)
  "cm":                               [variables.get_cm(mu_no) for mu_no in range(variables.n_compartments)],   # Cm parameter for every motor unit (capacitance of the cellular membrane)
  "timeStepWidth":                    variables.dt_0D,                    # time step width of the subcellular problem
  "endTime":                          variables.end_time,                 # end time, this is not relevant because it will be overridden by the splitting scheme
  "timeStepOutputInterval":           100,                                # how often the output timestep should be printed
  "solverName":                       "activationSolver",                 # reference to the solver used for the global linear system of the multidomain eq.
  "initialGuessNonzero":              True,                               # if the initial guess for the 3D system should be set as the solution of the previous timestep, this only makes sense for iterative solvers
  "inputIsGlobal":                    True,                               # if values and dofs correspond to the global numbering
  "showLinearSolverOutput":           False,                              # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
  "compartmentRelativeFactors":       variables.relative_factors.tolist(),          # list of lists of the factors for every dof, because "inputIsGlobal": True, this contains the global dofs
  "theta":                            0.5,                                # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
  "PotentialFlow": {
    "FiniteElementMethod" : {  
      "meshName":                     "3Dmesh",
      "solverName":                   "potentialFlowSolver",
      "prefactor":                    1.0,
      "dirichletBoundaryConditions":  variables.potential_flow_dirichlet_bc,
      "neumannBoundaryConditions":    [],
      "inputMeshIsGlobal":            True,
    },
  },
  "Activation": {
    "FiniteElementMethod" : {  
      "meshName":                     "3Dmesh",
      "solverName":                   "activationSolver",
      "prefactor":                    1.0,
      "inputMeshIsGlobal":            True,
      "dirichletBoundaryConditions":  {},
      "neumannBoundaryConditions":    [],
      "diffusionTensor": [[      # sigma_i           # fiber direction is (1,0,0)
        8.93, 0, 0,
        0, 0.0, 0,
        0, 0, 0.0
      ]], 
      "extracellularDiffusionTensor": [[      # sigma_e
        6.7, 0, 0,
        0, 6.7, 0,
        0, 0, 6.7,
      ]],
    },
  },
  "Fat": {
    "FiniteElementMethod" : {  
      "meshName":                     "3DFatMesh",
      "solverName":                   "activationSolver",
      "prefactor":                    0.4,
      "inputMeshIsGlobal":            True,
      "dirichletBoundaryConditions":  {},
      "neumannBoundaryConditions":    [],
    },
  },
  
  "OutputWriter" : [
    {"format": "Paraview", "outputInterval": (int)(1./variables.dt_0D*variables.output_timestep), "filename": "out/output", "binary": True, "fixedFormat": False, "combineFiles": True},
    #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
    #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
  ]
}
  
config = {
  "scenarioName":          variables.scenario_name,
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "meta": {                 # additional fields that will appear in the log
    "partitioning":         [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": variables.mappings_between_meshes,
  "Solvers": {
    "potentialFlowSolver": {
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e5,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFormat":         "default",
      "dumpFilename":       "",
    },
    "activationSolver": {
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e3,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFormat":         "matlab",
      "dumpFilename":       "",
    }
  },
  "StrangSplitting": {
    "timeStepWidth":          variables.dt_splitting,
    "logTimeStepWidthAsKey":  "dt_splitting",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 100,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": [0],          # CellML V_mk (0) <=> Multidomain V_mk^(i) (0)
    "connectedSlotsTerm2To1": [None, 0],    # Multidomain V_mk^(i+1) (1) -> CellML V_mk (0)

    "Term1": {      # CellML
      "MultipleInstances": {
        "nInstances": variables.n_compartments,  
        "instances": [        # settings for each motor unit, `i` is the index of the motor unit
        {
          "ranks": list(range(n_ranks)),
          "Heun" : {
            "timeStepWidth":                variables.dt_0D,  # 5e-5
            "logTimeStepWidthAsKey":        "dt_0D",
            "durationLogKey":               "duration_0D",
            "initialValues":                [],
            "timeStepOutputInterval":       1e4,
            "inputMeshIsGlobal":            True,
            "dirichletBoundaryConditions":  {},
            "nAdditionalFieldVariables":    0,
                
            "CellML" : {
              "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
              "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
              "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
              
              # optimization parameters
              "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
              "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
              "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
              "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
              
              # stimulation callbacks
              "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
              "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
              "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(compartment_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
              "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(compartment_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
              "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
              "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(compartment_no),# [ms] first time when to call setSpecificStates
              "additionalArgument":                     compartment_no,
              
              "mappings":                               variables.mappings,                             # mappings between parameters and intermediates/constants and between outputConnectorSlots and states, intermediates or parameters, they are defined in helper.py
              "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
              
              "meshName":                               "3Dmesh",
              "stimulationLogFilename":                 "out/stimulation.log",
            }
          }
        } for compartment_no in range(variables.n_compartments)]
      },
    },
    "Term2": {     # Diffusion, i.e. Multidomain
      "MultidomainSolver" : multidomain_solver,
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": (int)(1./variables.dt_splitting*variables.output_timestep), "filename": "out/surface", "binary": True, "fixedFormat": False, "combineFiles": True},
        ],
        "face": "1-",
        "MultidomainSolver" : multidomain_solver,
      }
    }
  }
}

print("Linear solver type: {}".format(config["Solvers"]["activationSolver"]["solverType"]))

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))

