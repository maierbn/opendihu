# Multiple 1D fibers (monodomain) with 3D intra-muscular EMG (static bidomain), biceps geometry
# to see all available arguments, execute: ./fibers_emg ../settings_fibers_emg.py -help
#
# if fiber_file=cuboid.bin, it uses a small cuboid test example
#
# You have to set n_subdomains such that it matches the number of processes, e.g. 2x2x1 = 4 processes.
# Decomposition is in x,y,z direction, the fibers are aligned with the z axis.
# E.g. --n_subdomains 2 2 1 which is 2x2x1 means no subdivision per fiber,
# --n_subdomains 8 8 4 means every fiber will be subdivided to 4 processes and all fibers will be computed by 8x8 processes.
#
# Example with 4 processes and end time 5, and otherwise default parameters:
#   mpirun -n 4 ./fibers_emg ../settings_fibers_emg.py --n_subdomains 2 2 1 --end_time=5.0
#
# Three files contribute to the settings:
# A lot of variables are set by the helper.py script, the variables and their values are defined in variables.py and this file
# creates the composite config that is needed by opendihu.
# You can provided parameter values in a custom_variables.py file in the variables subfolder of fibers_emg. (Instead of custom_variables.py you can choose any filename.)
# This custom variables file should be the next argument on the command line after settings_fibers_emg.py, e.g.:
#
#  ./fibers_emg ../settings_fibers_emg.py custom_variables.py --n_subdomains 1 1 1 --end_time=5.0
#  ./fibers_febio ../settings_fibers_emg.py febio.py

import sys, os
import timeit
import argparse
import importlib
import distutils.util

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
mbool = lambda x:bool(distutils.util.strtobool(x))   # function to parse bool arguments
parser = argparse.ArgumentParser(description='fibers_emg')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',            default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',             type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',                 type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',                 type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',                 type=int, default=variables.n_subdomains_z)
parser.add_argument('--diffusion_solver_type',               help='The solver for the diffusion.',                        default=variables.diffusion_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--diffusion_preconditioner_type',       help='The preconditioner for the diffusion.',                default=variables.diffusion_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--diffusion_solver_reltol',             help='Ralative tolerance for diffusion solver',              type=float, default=variables.diffusion_solver_reltol)
parser.add_argument('--diffusion_solver_maxit',              help='Maximum number of iterations for diffusion solver',    type=int, default=variables.diffusion_solver_maxit)
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).',  default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',           default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--potential_flow_solver_maxit',         help='Maximum number of iterations for potential flow solver', type=int, default=variables.potential_flow_solver_maxit)
parser.add_argument('--potential_flow_solver_reltol',        help='Relative tolerance for potential flow solver',         type=float, default=variables.potential_flow_solver_reltol)
parser.add_argument('--emg_solver_type',                     help='The solver for the static bidomain.',                  default=variables.emg_solver_type)
#parser.add_argument('--emg_solver_type',                    help='The solver for the static bidomain.',                  default=variables.emg_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--emg_preconditioner_type',             help='The preconditioner for the static bidomain.',          default=variables.emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_maxit',                    help='Maximum number of iterations for activation solver',   type=int, default=variables.emg_solver_maxit)
parser.add_argument('--emg_solver_reltol',                   help='Relative tolerance for activation solver',             type=float, default=variables.emg_solver_reltol)
parser.add_argument('--emg_solver_abstol',                   help='Absolute tolerance for activation solver',             type=float, default=variables.emg_solver_abstol)
parser.add_argument('--emg_initial_guess_nonzero',           help='If the initial guess for the emg linear system should be set to the previous solution.', default=variables.emg_initial_guess_nonzero, action='store_true')
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--python_output',                       help='Enable the python output writer.',                     default=variables.python_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--cellml_file',                         help='The filename of the file that contains the cellml model.', default=variables.cellml_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--stimulation_frequency',               help='Stimulations per ms. Each stimulation corresponds to one line in the firing_times_file.', default=variables.stimulation_frequency)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                             type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',                    type=float, default=variables.output_timestep)
parser.add_argument('--output_timestep_fibers',              help='The timestep for writing fiber outputs.',              type=float, default=variables.output_timestep_fibers)
parser.add_argument('--output_timestep_3D_emg',              help='The timestep for writing 3D emg outputs.',             type=float, default=variables.output_timestep_3D_emg)
parser.add_argument('--output_timestep_surface',             help='The timestep for writing 2D emg surface outputs.',     type=float, default=variables.output_timestep_surface)
parser.add_argument('--dt_0D',                               help='The timestep for the 0D model.',                       type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D',                               help='The timestep for the 1D model.',                       type=float, default=variables.dt_1D)
parser.add_argument('--dt_splitting',                        help='The timestep for the splitting.',                      type=float, default=variables.dt_splitting)
parser.add_argument('--dt_3D',                               help='The timestep for the 3D model, either bidomain or mechanics.', type=float, default=variables.dt_3D)
parser.add_argument('--optimization_type',                   help='The optimization_type in the cellml adapter.',         default=variables.optimization_type, choices=["vc", "simd", "openmp", "gpu"])
parser.add_argument('--approximate_exponential_function',    help='Approximate the exp function by a Taylor series',      type=mbool, default=variables.approximate_exponential_function)
parser.add_argument('--maximum_number_of_threads',           help='If optimization_type is "openmp", the max thread number, 0=all.', type=int, default=variables.maximum_number_of_threads)
parser.add_argument('--use_aovs_memory_layout',              help='If optimization_type is "vc", whether to use AoVS memory layout.', type=mbool, default=variables.use_aovs_memory_layout)
parser.add_argument('--disable_firing_output',               help='Disables the initial list of fiber firings.',          default=variables.disable_firing_output, action='store_true')
parser.add_argument('--enable_surface_emg',                  help='Enable the surface emg output writer.',                type=mbool, default=variables.enable_surface_emg)
parser.add_argument('--fast_monodomain_solver_optimizations',help='Enable the optimizations for fibers.',                 type=mbool, default=variables.fast_monodomain_solver_optimizations)
parser.add_argument('--enable_weak_scaling',                 help='Disable optimization for not stimulated fibers.',      default=False, action='store_true')
parser.add_argument('--v',                                   help='Enable full verbosity in c++ code')
parser.add_argument('-v',                                    help='Enable verbosity level in c++ code',                   action="store_true")
parser.add_argument('-vmodule',                              help='Enable verbosity level for given file in c++ code')
parser.add_argument('-pause',                                help='Stop at parallel debugging barrier',                   action="store_true")
parser.add_argument('--use_elasticity',                      help='Enable linear elasticity',                             action="store_true")
# parameter for the 3D mesh generation
parser.add_argument('--mesh3D_sampling_stride', nargs=3,     help='Stride to select the mesh points in x, y and z direction.', type=int, default=None)
parser.add_argument('--mesh3D_sampling_stride_x',            help='Stride to select the mesh points in x direction.',     type=int, default=variables.sampling_stride_x)
parser.add_argument('--mesh3D_sampling_stride_y',            help='Stride to select the mesh points in y direction.',     type=int, default=variables.sampling_stride_y)
parser.add_argument('--mesh3D_sampling_stride_z',            help='Stride to select the mesh points in z direction.',     type=int, default=variables.sampling_stride_z)

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
if variables.enable_weak_scaling:
  variables.fast_monodomain_solver_optimizations = False

# 3D mesh resolution
if variables.mesh3D_sampling_stride is not None:
    variables.mesh3D_sampling_stride_x = variables.mesh3D_sampling_stride[0]
    variables.mesh3D_sampling_stride_y = variables.mesh3D_sampling_stride[1]
    variables.mesh3D_sampling_stride_z = variables.mesh3D_sampling_stride[2]
variables.sampling_stride_x = variables.mesh3D_sampling_stride_x
variables.sampling_stride_y = variables.mesh3D_sampling_stride_y
variables.sampling_stride_z = variables.mesh3D_sampling_stride_z

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

if variables.use_elasticity:
  variables.emg_solver_type = "cg"
  
# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.1e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.1e}, potential_flow_solver_type: {}, approx. exp.: {}".format(variables.dt_1D, variables.potential_flow_solver_type, variables.approximate_exponential_function))
  print("dt_splitting:    {:0.1e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.dt_splitting, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.1e}, paraview_output: {}, optimization_type: {}{}".format(variables.dt_3D, variables.paraview_output, variables.optimization_type, " ({} threads)".format(variables.maximum_number_of_threads) if variables.optimization_type=="openmp" else " (AoVS)" if variables.optimization_type=="vc" and variables.use_aovs_memory_layout else " (SoVA)" if variables.optimization_type=="vc" and not variables.use_aovs_memory_layout else ""))
  print("output_timestep: {:0.1e}, surface: {:0.1e}, stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.output_timestep_surface, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("                          fast_monodomain_solver_optimizations: {}".format(variables.fast_monodomain_solver_optimizations))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  print("prefactor: sigma_eff/(Am*Cm) = {} = {} / ({}*{})".format(variables.Conductivity/(variables.Am*variables.Cm), variables.Conductivity, variables.Am, variables.Cm))
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# define the config dict
config = {
  "scenarioName":                   variables.scenario_name,    # scenario name which will appear in the log file
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/" + variables.scenario_name + "/mappings_between_meshes.txt",
  "meta": {                 # additional fields that will appear in the log
    "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {"MeshFiber_{}".format(i) : ["3Dmesh", "3Dmesh_quadratic"] for i in range(variables.n_fibers_total)},
  "Solvers": {
    "diffusionTermSolver": {# solver for the implicit timestepping scheme of the diffusion time step
      "relativeTolerance":  variables.diffusion_solver_reltol,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      variables.diffusion_solver_maxit,
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    },
    "potentialFlowSolver": {# solver for the initial potential flow, that is needed to estimate fiber directions for the bidomain equation
      "relativeTolerance":  variables.potential_flow_solver_reltol,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      variables.potential_flow_solver_maxit,
      "solverType":         variables.potential_flow_solver_type,
      "preconditionerType": variables.potential_flow_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
      "cycleType":          "cycleV",     # if the preconditionerType is "gamg", which cycle to use "cycleV" or "cycleW"
      "gamgType":           "agg",        # if the preconditionerType is "gamg", the type of the amg solver
      "nLevels":            25,           # if the preconditionerType is "gamg", the maximum number of levels
    },
    "muscularEMGSolver": {   # solver for the static Bidomain equation and the EMG
      "relativeTolerance":  variables.emg_solver_reltol,
      "absoluteTolerance":  variables.emg_solver_abstol,    
      "maxIterations":      variables.emg_solver_maxit,
      "solverType":         variables.emg_solver_type,
      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "linearElasticitySolver": {   # solver for linear elasticity
      "relativeTolerance":  1e-1,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    ,
      "maxIterations":      1e4,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    }, 
  },
  "Coupling": {
    "timeStepWidth":          variables.dt_3D,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_3D",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": {1:0} if variables.use_elasticity else {0:0},  # elasticity: transfer gamma to elasticity, fibers_emg: transfer Vm to StaticBidomainSolver
    "connectedSlotsTerm2To1": [None]   if variables.use_elasticity else [None],   # elasticity: only transfer back geometry (this happens automatically),   fibers_emg: transfer nothing back
    "Term1": {        # monodomain, fibers
      "MultipleInstances": {
        "logKey":                     "duration_subdomains_xy",
        "ranksAllComputedInstances":  list(range(n_ranks)),
        "nInstances":                 variables.n_subdomains_xy,
        "instances": 
        [{
          "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
          "StrangSplitting": {
            #"numberTimeSteps": 1,
            "timeStepWidth":          variables.dt_splitting,  # 1e-1
            "logTimeStepWidthAsKey":  "dt_splitting",
            "durationLogKey":         "duration_monodomain",
            "timeStepOutputInterval": 100,
            "endTime":                variables.dt_splitting,
            "connectedSlotsTerm1To2": [0,1] if variables.use_elasticity else [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion), for elasticity also transfer gamma
            "connectedSlotsTerm2To1": [0,1] if variables.use_elasticity else [0],   # transfer the same back, this avoids data copy

            "Term1": {      # CellML, i.e. reaction term of Monodomain equation
              "MultipleInstances": {
                "logKey":             "duration_subdomains_z",
                "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                "instances": 
                [{
                  "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "Heun" : {
                    "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                    "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                    "durationLogKey":               "duration_0D",                           # log key of duration for this solver
                    "timeStepOutputInterval":       1e4,                                     # how often to print the current timestep
                    "initialValues":                [],                                      # no initial values are specified
                    "dirichletBoundaryConditions":  {},                                      # no Dirichlet boundary conditions are specified
                    "dirichletOutputFilename":      None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                    
                    "inputMeshIsGlobal":            True,                                    # the boundary conditions and initial values would be given as global numbers
                    "checkForNanInf":               False,                                   # abort execution if the solution contains nan or inf values
                    "nAdditionalFieldVariables":    0,                                       # number of additional field variables
                    "additionalSlotNames":          [],                                      # names for the additional slots
                      
                    "CellML" : {
                      "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                      #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                      "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                      "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                      
                      # optimization parameters
                      "optimizationType":                       variables.optimization_type,                    # "vc", "simd", "openmp" type of generated optimizated source file
                      "approximateExponentialFunction":         variables.approximate_exponential_function,     # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                      "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                      "maximumNumberOfThreads":                 variables.maximum_number_of_threads,            # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                      "useAoVSMemoryLayout":                    variables.use_aovs_memory_layout,               # if optimizationType is "vc", whether to use the Array-of-Vectorized-Struct (AoVS) memory layout instead of the Struct-of-Vectorized-Array (SoVA) memory layout. Setting to True is faster.
                      
                      # stimulation callbacks
                      #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                      #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                      #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                      #"setSpecificStatesCallInterval":          2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                      "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                      "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                      "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                      "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                      "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                      
                      # parameters to the cellml model
                      "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                      "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                      
                      "meshName":                               "MeshFiber_{}".format(fiber_no),                # reference to the fiber mesh
                      "stimulationLogFilename":                 "out/" + variables.scenario_name + "/stimulation.log",                          # a file that will contain the times of stimulations
                    },      
                    "OutputWriter" : [
                      {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
                    ] if variables.states_output else []
                    
                  },
                } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                    for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                      for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                        for motor_unit_no in [get_motor_unit_no(fiber_no)]],
              }
            },
            "Term2": {     # Diffusion
              "MultipleInstances": {
                "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                "instances": 
                [{
                  "ranks":                         list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "ImplicitEuler" : {       # include both CrankNicolson and ImplicitEuler in the settings such that both variants in C++ file are possible
                    "initialValues":               [],                                      # initial values to be set in the solution vector prior to the first timestep
                    #"numberTimeSteps":            1,
                    "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                    "timeStepWidthRelativeTolerance": 1e-10,                                # tolerance for the time step width, when to rebuild the system matrix
                    "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                    "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                    "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep
                    "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                    "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                    "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                    "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                    "checkForNanInf":              False,                                   # if the solution should be checked for NaN and inf values, this requires a lot of runtimes
                    "nAdditionalFieldVariables":   1 if variables.use_elasticity else 2,    # number of additional field variables that should be added and potentially written to output files, these field variables can be used for receiving data from other solvers
                    "additionalSlotNames":         [],                                      # slot names for the additional field variables
                    "FiniteElementMethod" : {
                      "inputMeshIsGlobal":         True,
                      "meshName":                  "MeshFiber_{}".format(fiber_no),
                      "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                      "solverName":                "diffusionTermSolver",
                      "slotName":                  "",
                    },
                    "OutputWriter" : [
                      #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                      #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                      #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                      #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                    ]
                  },
                  "CrankNicolson": {      # include both CrankNicolson and ImplicitEuler in the settings such that both variants in C++ file are possible
                    "initialValues":               [],                                      # initial values to be set in the solution vector prior to the first timestep
                    #"numberTimeSteps":            1,
                    "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                    "timeStepWidthRelativeTolerance": 1e-10,                                # tolerance for the time step width, when to rebuild the system matrix
                    "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                    "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                    "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep to console
                    "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                    "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                    "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                    "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                    "checkForNanInf":              False,                                   # if the solution should be checked for NaN and inf values, this requires a lot of runtimes
                    "nAdditionalFieldVariables":   1 if variables.use_elasticity else 0,    # number of additional field variables that should be added and potentially written to output files, these field variables can be used for receiving data from other solvers
                    "additionalSlotNames":         [],                                      # slot names for the additional field variables
                    "FiniteElementMethod" : { 
                      "inputMeshIsGlobal":         True,
                      "meshName":                  "MeshFiber_{}".format(fiber_no),
                      "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                      "solverName":                "diffusionTermSolver",
                      "slotName":                  "", 
                    },  
                    "OutputWriter" : [ 
                    ]   
                  },   
                } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                    for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                      for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                        for motor_unit_no in [get_motor_unit_no(fiber_no)]],
                "OutputWriter" : variables.output_writer_fibers,
              },
            },
          }
        } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
        for subdomain_coordinate_y in range(variables.n_subdomains_y)
            for subdomain_coordinate_x in range(variables.n_subdomains_x)]
      },
      "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
      "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
      "onlyComputeIfHasBeenStimulated": variables.fast_monodomain_solver_optimizations,                          # only compute fibers after they have been stimulated for the first time
      "disableComputationWhenStatesAreCloseToEquilibrium": variables.fast_monodomain_solver_optimizations,       # optimization where states that are close to their equilibrium will not be computed again      
      "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set
      "neuromuscularJunctionRelativeSize": 0.1,                        # range where the neuromuscular junction is located around the center, relative to fiber length. The actual position is draws randomly from the interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at the center, 0.1 means located approximately at the center, but it can vary 10% in total between all fibers.
      "generateGPUSource":        True,                                # (set to True) only effective if optimizationType=="gpu", whether the source code for the GPU should be generated. If False, an existing source code file (which has to have the correct name) is used and compiled, i.e. the code generator is bypassed. This is useful for debugging, such that you can adjust the source code yourself. (You can also add "-g -save-temps " to compilerFlags under CellMLAdapter)
      "useSinglePrecision":       False,                               # only effective if optimizationType=="gpu", whether single precision computation should be used on the GPU. Some GPUs have poor double precision performance. Note, this drastically increases the error and, in consequence, the timestep widths should be reduced.
      #"preCompileCommand":        "bash -c 'module load argon-tesla/gcc/11-20210110-openmp; module list; gcc --version",     # only effective if optimizationType=="gpu", system command to be executed right before the compilation
      #"postCompileCommand":       "'",   # only effective if optimizationType=="gpu", system command to be executed right after the compilation
    },
    "Term2": {        # Bidomain, EMG
      "StaticBidomainSolver": {       # version for fibers_emg
        "timeStepWidth":          variables.dt_3D,
        "timeStepOutputInterval": 50,
        "durationLogKey":         "duration_bidomain",
        "solverName":             "muscularEMGSolver",
        "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
        "slotNames":              [],
        
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
            "solverName":         "muscularEMGSolver",
            "prefactor":          1.0,
            "inputMeshIsGlobal":  True,
            "dirichletBoundaryConditions": {},
            "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            "neumannBoundaryConditions":   [],
            "slotName":           "",
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
      },
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_surface), "filename": "out/" + variables.scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental",},
        ] if variables.enable_surface_emg else [],
        "face":                     ["0+"],              # which faces of the 3D mesh should be written into the 2D mesh
        "samplingPoints":           None,                # the electrode positions, they are created in the helper.py script
        "updatePointPositions":     False,               # the electrode points should be initialize in every timestep (set to False for the static case). This makes a difference if the muscle contracts, then True=fixed electrodes, False=electrodes moving with muscle.
        "filename":                 "out/{}/electrodes.csv".format(variables.scenario_name),
        "enableCsvFile":            False,               # if the values at the sampling points should be written to csv files
        "enableVtpFile":            False,               # if the values at the sampling points should be written to vtp files
        "enableGeometryInCsvFile":  False,               # if the csv output file should contain geometry of the electrodes in every time step. This increases the file size and only makes sense if the geometry changed throughout time, i.e. when computing with contraction
        "enableGeometryFiles":      False,               # if there should be extra files of the locations of the electrodes on every rank
        "xiTolerance":              0.3,                 # tolerance for element-local coordinates xi, for finding electrode positions inside the elements. Increase or decrease this numbers if not all electrode points are found.
        
        "StaticBidomainSolver": {
          "timeStepWidth":          variables.dt_3D,
          "timeStepOutputInterval": 50,
          "durationLogKey":         "duration_bidomain",
          "solverName":             "muscularEMGSolver",
          "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
          "enableJacobianConditionNumber": False,        # if set to true, estimate the condition number of the jacobian of the element-coordinate-to-world-frame mapping in every element and output it in the output writer
          "slotNames":              [],

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
              "solverName":         "muscularEMGSolver",
              "prefactor":          1.0,
              "inputMeshIsGlobal":  True,
              "dirichletBoundaryConditions": {},
              "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "neumannBoundaryConditions":   [],
              "slotName":           "",
              "diffusionTensor": [[      # sigma_i, fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element
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
      },
      "QuasiStaticLinearElasticitySolver": {
        "PotentialFlow": {        # potential flow for fiber directions in the 3D mesh
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
        "FiniteElementMethod" : {   # linear elasticity finite element method
          "meshName":             "3Dmesh",
          "solverName":           "linearElasticitySolver",
          "prefactor":            1.0,
          "inputMeshIsGlobal":    True,
          "dirichletBoundaryConditions": variables.use_elasticity_dirichlet_bc,
          "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
          "neumannBoundaryConditions":   variables.use_elasticity_neumann_bc,
          "bulkModulus":          40e3, #40e3 # https://www.researchgate.net/publication/230248067_Bulk_Modulus
          "shearModulus":         39e3, #39e3 # https://onlinelibrary.wiley.com/doi/full/10.1002/mus.24104
          "slotName":             "",
        },
        "maximumActiveStress":      1.0,
        "strainScalingCurveWidth":  1.0,
        "scalingFactor":            1e4,   #1e4
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/deformation", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
          #{"format": "PythonFile", "filename": "out/deformation", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
        ]
      },
      "QuasiStaticNonlinearElasticitySolverFebio": {
        "durationLogKey": "febio",
        "meshName":       "3Dmesh",
        "activationFactor": 1e-5,
        "force":            100,     # force on top of muscle
        "OutputWriter" : variables.output_writer_elasticity,
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
