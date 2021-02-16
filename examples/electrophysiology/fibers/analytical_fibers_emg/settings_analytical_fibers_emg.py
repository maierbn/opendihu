# Multiple 1D fibers with analytical Vm computation with 3D EMG (static bidomain), biceps geometry
# to see all available arguments, execute: ./analytical_fibers_emg ../settings_fibers_emg.py -help
#

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

# define command line arguments
parser = argparse.ArgumentParser(description='fibers_emg')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',            default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',             type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',                 type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',                 type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',                 type=int, default=variables.n_subdomains_z)
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).',  default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',           default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--potential_flow_solver_maxit',         help='Maximum number of iterations for potential flow solver', type=int, default=variables.potential_flow_solver_maxit)
parser.add_argument('--potential_flow_solver_reltol',        help='Relative tolerance for potential flow solver',         type=float, default=variables.potential_flow_solver_reltol)
parser.add_argument('--emg_solver_type',                     help='The solver for the static bidomain.',                  default=variables.emg_solver_type)
#parser.add_argument('--emg_solver_type',                    help='The solver for the static bidomain.',                  default=variables.emg_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--emg_preconditioner_type',             help='The preconditioner for the static bidomain.',          default=variables.emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_maxit',                    help='Maximum number of iterations for activation solver',   type=int, default=variables.emg_solver_maxit)
parser.add_argument('--emg_solver_reltol',                   help='Ralative tolerance for activation solver',             type=float, default=variables.emg_solver_reltol)
parser.add_argument('--emg_initial_guess_nonzero',           help='If the initial guess for the emg linear system should be set to the previous solution.', default=variables.emg_initial_guess_nonzero, action='store_true')
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',                   default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',              default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--stimulation_frequency',               help='Stimulations per ms. Each stimulation corresponds to one line in the firing_times_file.', default=variables.stimulation_frequency)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                             type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',                    type=float, default=variables.output_timestep)
parser.add_argument('--dt_3D',                               help='The timestep for the 3D model, either bidomain or mechanics.', type=float, default=variables.dt_3D)
parser.add_argument('--disable_firing_output',               help='Disables the initial list of fiber firings.',          default=variables.disable_firing_output, action='store_true')
parser.add_argument('--enable_surface_emg',                  help='Enable the surface emg output writer.',                default=variables.enable_surface_emg, action='store_true')
parser.add_argument('--v',                                   help='Enable full verbosity in c++ code')
parser.add_argument('-v',                                    help='Enable verbosity level in c++ code',                   action="store_true")
parser.add_argument('-vmodule',                              help='Enable verbosity level for given file in c++ code')
parser.add_argument('-on_error_attach_debugger',             help='Enable verbosity level for given file in c++ code')
parser.add_argument('-pause',                                help='Stop at parallel debugging barrier',                   action="store_true")

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
    for j in range(1,n_ranks+1):
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
  print("potential_flow_solver_type: {}, emg_solver_type: {}, emg_initial_guess_nonzero: {}".format(variables.potential_flow_solver_type, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

if variables.custom_meshes:
  variables.meshes = variables.custom_meshes
  variables.n_fibers_x = variables.n_custom_fibers_x
  variables.n_fibers_y = variables.n_custom_fibers_y
  variables.n_fibers_per_subdomain_x = variables.n_fibers_x
  variables.n_fibers_per_subdomain_y = variables.n_fibers_y
  variables.n_fibers_total = variables.n_custom_fibers
  
  print("set n_fibers_total: {}".format(variables.n_fibers_total))

#print("variables.meshes: {}".format(variables.meshes))

# define the config dict
config = {
  "scenarioName":                   variables.scenario_name,    # scenario name which will appear in the log file
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",     # output file that contains a log about creation of mappings between meshes
  "meta": {                 # additional fields that will appear in the log
    "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes": variables.meshes,
  "MappingsBetweenMeshes": {
    "MeshFiber_{}".format(i): {
      "name": "3Dmesh", 
      "xiTolerance": 0.1, 
      "defaultValue": -90,
      "enableWarnings": True,
      "compositeUseOnlyInitializedMappings": False,
      "fixUnmappedDofs": True, 
    } for i in range(variables.n_fibers_total)
  },
  "Solvers": {
    "potentialFlowSolver": {# solver for the initial potential flow, that is needed to estimate fiber directions for the bidomain equation
      "relativeTolerance":  variables.potential_flow_solver_reltol,
      "absoluteTolerance":  1e-4,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      variables.potential_flow_solver_maxit,
      "solverType":         variables.potential_flow_solver_type,
      "preconditionerType": variables.potential_flow_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "activationSolver": {   # solver for the static Bidomain equation and the EMG
      "relativeTolerance":  variables.emg_solver_reltol,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      variables.emg_solver_maxit,
      "solverType":         variables.emg_solver_type,
      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "ascii",
    },
  },
  "Coupling": {
    "timeStepWidth":          variables.dt_3D,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_3D",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 10,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": None,
    "connectedSlotsTerm2To1": None,
    "Term1": {        # analytical fibers
      "MultipleInstances": {
        "logKey":                     "duration_fibers",
        "ranksAllComputedInstances":  list(range(n_ranks)),
        "nInstances":                 variables.n_fibers_total,
        "instances":
        [{
          "ranks":                    list(range(subdomain_no,n_ranks,variables.n_subdomains_xy)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
          "PrescribedValues": {
            "meshName":               "MeshFiber_{}".format(fiber_no),               # reference to the fiber mesh
            "numberTimeSteps":        1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
            "timeStepOutputInterval": 20,            # if the time step should be written to console, a value > 10 produces no output
            "slotNames":              ["vm"],        # names of the connector slots
            
            # a list of field variables that will get values assigned in every timestep, by the provided callback function
            "fieldVariables1": [
              {"name": "Vm",     "callback": set_vm_values},
            ],
            "fieldVariables2":     [],
            "additionalArgument":  fiber_no,         # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
            
          },      
        } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
        for subdomain_coordinate_y in range(variables.n_subdomains_y) \
          for subdomain_coordinate_x in range(variables.n_subdomains_x) \
            for subdomain_no in [subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x] \
              for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                  for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                    for motor_unit_no in [get_motor_unit_no(fiber_no)]],
                    
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_fibers), "filename": "out/" + variables.scenario_name + "/prescribed_fiber", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
        ]
      }
    },
    "Term2": {        # Bidomain, EMG
      "StaticBidomainSolver": {       # version for fibers_emg
        "timeStepWidth":          variables.dt_3D,
        "timeStepOutputInterval": 50,
        "durationLogKey":         "duration_bidomain",
        "solverName":             "activationSolver",
        "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
        "slotNames":              ["vm"],
        
        "PotentialFlow": {
          "FiniteElementMethod" : {
            "meshName":           "3Dmesh",
            "solverName":         "potentialFlowSolver",
            "prefactor":          1.0,
            "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
            "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            "neumannBoundaryConditions":   [],
            "inputMeshIsGlobal":  True,
            "slotName":           "",
          },
        },
        "Activation": {
          "FiniteElementMethod" : {
            "meshName":           "3Dmesh",
            "solverName":         "activationSolver",
            "prefactor":          1.0,
            "inputMeshIsGlobal":  True,
            "dirichletBoundaryConditions": {},
            "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            "neumannBoundaryConditions":   [],
            "slotName":           "vm",
            "diffusionTensor": [[      # sigma_i,  fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element
              8.93, 0, 0,
              0, 0.893, 0,
              0, 0, 0.893
            ]],
            "extracellularDiffusionTensor": [[      # sigma_e, one list item = same tensor for all elements, multiple list items = a different tensor for each element
              6.7, 0, 0,
              0, 6.7, 0,
              0, 0, 6.7,
            ]],
          },
        },
        "OutputWriter" : variables.output_writer_emg,
      }
    }
  }
}

if False:
  instances = config["Coupling"]["Term1"]["MultipleInstances"]["instances"]
  for i,instance in enumerate(instances):
    if instance is not None:
      print("{}: instance {}, ranks: {}".format(rank_no,i,instance["ranks"]))

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
