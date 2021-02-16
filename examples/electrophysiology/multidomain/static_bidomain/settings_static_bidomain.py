# Multiple 1D fibers (monodomain) with 3D intra-muscular EMG (static bidomain) and 3D fat layer (anisotropic diffusion), on biceps geometry

import sys, os
import timeit
import argparse
import importlib

# parse rank arguments
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
    print("Warning: There is no variables file, e.g:\n ./prescribed_static_bidomain_solver ../settings_static_bidomain.py normal.py\n")
  exit(0)

# -------------- begin user parameters ----------------

# -------------- end user parameters ----------------

# define command line arguments
parser = argparse.ArgumentParser(description='prescribed_static_bidomain_solver')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).', default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',  default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_type',                     help='The solver for the static bidomain.',         default=variables.emg_solver_type)
parser.add_argument('--emg_preconditioner_type',             help='The preconditioner for the static bidomain.', default=variables.emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_initial_guess_nonzero',           help='If the initial guess for the emg linear system should be set to the previous solution.', default=variables.emg_initial_guess_nonzero, action='store_true')
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',           type=float, default=variables.output_timestep)
parser.add_argument('--dt_coupling',                         help='The timestep for the coupling.',             type=float, default=variables.dt_coupling)
parser.add_argument('--v',                                   help='Enable full verbosity in c++ code')
parser.add_argument('-v',                                    help='Enable verbosity level in c++ code', action="store_true")
parser.add_argument('-vmodule',                              help='Enable verbosity level for given file in c++ code')
parser.add_argument('-pause',                                help='Stop at parallel debugging barrier', action="store_true")

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

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
  for i in range(1,n_ranks+1):
    for j in [1]:
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = (int)(n_ranks / (i*j))
        performance = (k-optimal_value)**2 + (j-optimal_value)**2 + 1.1*(i-optimal_value)**2
        possible_partitionings.append([i,j,k,performance])
        
  # if no possible partitioning was found
  if len(possible_partitionings) == 0:
    if rank_no == 0:
      print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {} and no automatic partitioning could be done.\n\n\033[0m".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
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
  print("dt_coupling:     {:0.0e}".format(variables.dt_coupling))
  print("potential_flow_solver_type: {}".format(variables.potential_flow_solver_type))
  print("emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

config_static_bidomain_solver = {             # solves Bidomain equation: K(sigma_i) Vm + K(sigma_i+sigma_e) phi_e = 0   => K(sigma_i+sigma_e) phi_e = -K(sigma_i) Vm
  "numberTimeSteps":        1,
  "timeStepOutputInterval": 50,
  "durationLogKey":         "duration_bidomain",
  "solverName":             "muscularEMGSolver",
  "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
  "slotNames":              ["vm"],
  
  "PotentialFlow": {
    "FiniteElementMethod" : {  
      "meshName":           ["3Dmesh","3DFatMesh"],
      "solverName":         "potentialFlowSolver",
      "prefactor":          1.0,
      "slotName":           "",
      "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
      "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "neumannBoundaryConditions":   [],
      "inputMeshIsGlobal":  True,
    },
  },
  "Activation": {
    "FiniteElementMethod" : {  
      "meshName":           ["3Dmesh","3DFatMesh"],       # composite mesh that consists of the muscle value mesh and the fat layer mesh
      "solverName":         "muscularEMGSolver",
      "prefactor":          1.0,
      "inputMeshIsGlobal":  True,
      "dirichletBoundaryConditions": {},
      "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "neumannBoundaryConditions":   [],
      "slotName":           "",
      
      # ∇•(sigma_i+sigma_e)∇phi_e = -∇•(sigma_i)∇Vm
      "diffusionTensor": [
        [                
          8.93, 0, 0,     # sigma_i, for muscle volume                # fiber direction is (1,0,0)
          0, 0.893, 0,
          0, 0, 0.893
        ],[ 
          0, 0, 0,        # no sigma_i for fat mesh!
          0, 0, 0,
          0, 0, 0
        ]
      ],
      "extracellularDiffusionTensor": [      # sigma_e
        [
          6.7, 0, 0,      # sigma_e, conductivity in extra-cellular space
          0, 6.7, 0,
          0, 0, 6.7,
        ],[
          0.4, 0, 0,      # sigma, conductivity in fat layer
          0, 0.4, 0,
          0, 0, 0.4,
        ],
      ]
    },
  },
  "OutputWriter" : variables.output_writer_emg,
}

# define the config dict
config = {
  "scenarioName":          variables.scenario_name,
  "mappingsBetweenMeshesLogFile":  "out/mappings_between_meshes.txt",
  "logFormat":             "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "meta": {                 # additional fields that will appear in the log
    "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": variables.mappings_between_meshes,
  "Solvers": {
    "potentialFlowSolver": {# solver for the initial potential flow, that is needed to estimate fiber directions for the bidomain equation
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      1e4,
      "solverType":         variables.potential_flow_solver_type,
      "preconditionerType": variables.potential_flow_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "muscularEMGSolver": {   # solver for the static Bidomain equation and the EMG
      "relativeTolerance":  1e-15,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      1e4,
      "solverType":         variables.emg_solver_type,
      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
  },
  "StaticBidomainSolver": config_static_bidomain_solver,        # this is for the executable `static_bidomain_solver`
  "Coupling": {                                                 # this is for the executable `prescribed_static_bidomain_solver`
    "timeStepWidth":          variables.dt_coupling,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_coupling",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": {0:0},
    "connectedSlotsTerm2To1": [],
    "Term1": {        # monodomain, fibers
      "PrescribedValues": {
        "meshName":                               "3Dmesh",                                       # use the linear mesh, it was partitioned by the helper.py script which called opendihu/scripts/create_partitioned_meshes_for_settings.py
        "additionalArgument":                     None,
        "slotNames":                              [],
        "numberTimeSteps":                        1,
        "timeStepOutputInterval":                 1,
        "fieldVariables1": [
          {"name": "Vm", "callback": variables.set_artifical_activation_values},
        ],
        "fieldVariables2":                        [],
      },
    },
    "Term2": {
      "StaticBidomainSolver": config_static_bidomain_solver,
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
