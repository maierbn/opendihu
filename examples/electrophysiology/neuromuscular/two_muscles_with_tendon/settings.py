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
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there  

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
parser.add_argument('--cellml_file',                         help='The filename of the file that contains the cellml model.', default=variables.cellml_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--stimulation_frequency',               help='Stimulations per ms. Each stimulation corresponds to one line in the firing_times_file.', default=variables.stimulation_frequency)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                             type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',                    type=float, default=variables.output_timestep)
parser.add_argument('--output_timestep_fibers',              help='The timestep for writing fiber outputs.',              type=float, default=variables.output_timestep_fibers)
parser.add_argument('--output_timestep_3D_emg',              help='The timestep for writing 3D emg outputs.',             type=float, default=variables.output_timestep_emg)
parser.add_argument('--output_timestep_surface',             help='The timestep for writing 2D emg surface outputs.',     type=float, default=variables.output_timestep_surface)
parser.add_argument('--dt_0D',                               help='The timestep for the 0D model.',                       type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D',                               help='The timestep for the 1D model.',                       type=float, default=variables.dt_1D)
parser.add_argument('--dt_splitting',                        help='The timestep for the splitting.',                      type=float, default=variables.dt_splitting_0D1D)
parser.add_argument('--dt_bidomain',                         help='The timestep for the 3D bidomain model.', type=float, default=variables.dt_bidomain)
parser.add_argument('--dt_elasticity',                       help='The timestep for the elasticity model.', type=float, default=variables.dt_elasticity)
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
        k = int(n_ranks / (i*j))
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
  print("dt_0D:           {:0.1e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.1e}, potential_flow_solver_type: {}, approx. exp.: {}".format(variables.dt_1D, variables.potential_flow_solver_type, variables.approximate_exponential_function))
  print("dt_splitting:    {:0.1e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.dt_splitting_0D1D, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.1e}, paraview_output: {}, optimization_type: {}{}".format(variables.dt_bidomain, variables.paraview_output, variables.optimization_type, " ({} threads)".format(variables.maximum_number_of_threads) if variables.optimization_type=="openmp" else " (AoVS)" if variables.optimization_type=="vc" and variables.use_aovs_memory_layout else " (SoVA)" if variables.optimization_type=="vc" and not variables.use_aovs_memory_layout else ""))
  print("output_timestep: {:0.1e}, surface: {:0.1e}, stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.output_timestep_surface, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("                          fast_monodomain_solver_optimizations: {}".format(variables.fast_monodomain_solver_optimizations))
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


# meshes

# add neuron meshes
neuron_meshes = {
  "motoneuronMesh": {
    "nElements" :         variables.n_motoneurons-1 if n_ranks == 1 else variables.n_motoneurons*n_ranks,  # the last dof is empty in parallel
    "physicalExtent":     1, # has no special meaning. only to seperate data points in paraview
    "physicalOffset":     0,
    "logKey":             "motoneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleMesh": {
    "nElements" :         variables.n_muscle_spindles-1 if n_ranks == 1 else variables.n_muscle_spindles*n_ranks,
    "physicalExtent":     1, # has no special meaning. only to seperate data points in paraview
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "golgiTendonOrganMesh": {
    "nElements" :         variables.n_golgi_tendon_organs-1 if n_ranks == 1 else variables.n_golgi_tendon_organs*n_ranks,
    "physicalExtent":     1, # has no special meaning. only to seperate data points in paraview
    "physicalOffset":     0,
    "logKey":             "golgi_tendon_organ",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  # no `nodePositions` fields as the nodes are created internally
  "muscle1Mesh": {
    "nElements" :         variables.n_elements_muscle1,
    "physicalExtent":     variables.muscle1_extent,
    "physicalOffset":     [0,0,0],
    "logKey":             "muscle1",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscle2Mesh": {
    "nElements" :         variables.n_elements_muscle2,
    "physicalExtent":     variables.muscle2_extent,
    "physicalOffset":     [0.0, 0.0, variables.muscle1_extent[2] + variables.tendon_length],
    "logKey":             "muscle2",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  # needed for mechanics solver
  "muscle1Mesh_quadratic": {
    "nElements" :         [elems // 2 for elems in variables.n_elements_muscle1],
    "physicalExtent":     variables.muscle1_extent,
    "physicalOffset":     [0,0,0],
    "logKey":             "muscle1_quadratic",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks,
  },
  "muscle2Mesh_quadratic": {
    "nElements" :         [elems // 2 for elems in variables.n_elements_muscle2],
    "physicalExtent":     variables.muscle2_extent,
    "physicalOffset":     [0.0, 0.0, variables.muscle1_extent[2] + variables.tendon_length],
    "logKey":             "muscle2_quadratic",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks,
  },
}
variables.meshes.update(neuron_meshes)
variables.meshes.update(fiber_meshes)


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
  # meshes know their coordinates and the mapping happens automatically. We could also define other parameters such as a mapping tolerance
  "MappingsBetweenMeshes": {"muscle{}_fiber{}".format(m,f) : ["muscle{}Mesh", "muscle{}Mesh_quadratic"] for f in range(variables.n_fibers_total) for m in [1,2]},
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
    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":   variables.linear_relative_tolerance,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   variables.linear_absolute_tolerance,           # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          variables.elasticity_solver_type,            # type of the linear solver
      "preconditionerType":  variables.elasticity_preconditioner_type,    # type of the preconditioner
      "maxIterations":       1e4,                                         # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                                  # maximum number of function iterations
      "snesMaxIterations":   variables.snes_max_iterations,               # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": variables.snes_relative_tolerance,         # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": variables.snes_absolute_tolerance,         # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",                                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": variables.snes_rebuild_jacobian_frequency,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "hypreOptions":        "",                                          # additional options for the hypre solvers could be given here
      "dumpFilename":        "",                                          # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",                                    # default, ascii, matlab
    }
  },


  # connections of the slots, identified by slot name
  "connectedSlots": [
    # global slots only support named slots (connectedSlotsTerm1To2 also allows indices)

    # use global slot, because automatic connection of "Razumova/activestress" does not work for some reason
    ("m1gout", "m1g_in"),
    ("m2gout", "m2g_in"),

    # lambda and derived values (by MapDofs) -> input of muscle splindel simulation
    ("m1ms0",    "ms_in0"),
    ("m1ms1",    "ms_in1"),
    ("m1ms2",    "ms_in2"),
    ("m1ms3",    "ms_in3"),
    ("m1ms4",    "ms_in4"),
  ],


  "Coupling": {
    "timeStepWidth":          variables.dt_elasticity,
    "logTimeStepWidthAsKey":  "dt_elasticity",
    "durationLogKey":         "duration_coupling",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    # "connectedSlotsTerm1To2": {"m1gout":"m1g_in", "m2gout":"m2g_in"}, # only use numbers here!
    # "connectedSlotsTerm1To2": {11:2,19:9},       # same effect as global slot definition
    "connectedSlotsTerm2To1": None,       # transfer nothing back

    "Term1": {
      "MultipleCoupling": {
       "timeStepWidth":          variables.end_time,
       "logTimeStepWidthAsKey":  "dt_multiple_coupling",
       "durationLogKey":         "duration_multiple_coupling",
       "timeStepOutputInterval": 1,
       "endTime":                variables.end_time,
       "deferInitialization":    True,       # initialize nested solvers only right before computation of first timestep
       "connectedSlotsTerm1To2": None,       # connect lambda to slot 0 and gamma to slot 2
       "connectedSlotsTerm2To1": None,       # transfer nothing back

        "Term1": {
          # spindles

          # mapping muscle spindles output -> motor neuron input
          "MapDofs": {
            "description":                "muscle_spindles_to_motoneurons",   # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                                  # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["mn_in"],
            "meshName":                   "motoneuronMesh",                   # the mesh on which the additional field variables will be defined
            "beforeComputation": None,
            "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
              {
                "fromConnectorSlot":                "ms_out",
                "toConnectorSlots":                 ["mn_in"],
                "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        list(range(variables.n_muscle_spindles)),   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleMesh"
                "outputDofs":                       [list(range(variables.n_motoneurons))],   # [0,1,...,n_motor_neurons], this is for mesh "motoneuronMesh"
                "callback":                         variables.callback_muscle_spindles_to_motoneurons,
                #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
            ],

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
                  {"format": "Paraview",   "outputInterval": int(2./variables.dt_muscle_spindles*variables.output_timestep_spindles), "filename": "out/" + variables.scenario_name + "/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                  {"format": "PythonFile", "outputInterval": int(2./variables.dt_muscle_spindles*variables.output_timestep_spindles), "filename": "out/" + variables.scenario_name + "/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                ]
              }
            }
          }
        },
        "Term2": {
          # Golgi tendon organ model solver
          "Heun": {
            "description":                  "Golgi tendon organs",
            "timeStepWidth":                variables.dt_golgi_tendon_organs,
            "logTimeStepWidthAsKey":        "dt_golgi_tendon_organs",
            "durationLogKey":               "duration_golgi_tendon_organ",
            "initialValues":                [],
            "timeStepOutputInterval":       500,
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
                {"format": "Paraview",   "outputInterval": int(2./variables.dt_golgi_tendon_organs*variables.output_timestep_golgi_tendon_organs), "filename": "out/" + variables.scenario_name + "/golgi_tendon_organs", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                {"format": "PythonFile", "outputInterval": int(2./variables.dt_golgi_tendon_organs*variables.output_timestep_golgi_tendon_organs), "filename": "out/" + variables.scenario_name + "/golgi_tendon_organs", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
              ]
            }
          }
        },
        "Term3": {
          # motor neurons
          "Heun" : {
            "description":                  "motoneurons",
            "timeStepWidth":                variables.dt_motoneuron,
            "logTimeStepWidthAsKey":        "dt_motoneuron",
            "durationLogKey":               "duration_motoneuron",
            "initialValues":                [],
            "timeStepOutputInterval":       500,
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
                {"format": "Paraview",   "outputInterval": int(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "filename": "out/" + variables.scenario_name + "/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                {"format": "PythonFile", "outputInterval": int(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "filename": "out/" + variables.scenario_name + "/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
              ]
            }
          }
        },
        "Term4": {
          # muscle1: bidoamin + 1D monodomain + 0D

          # TODO in other examples MapDofs also wraps the coupled mechanics solver
          # map from λ in the 3D mesh to muscle spindles input
          "MapDofs": {
            "description":                "muscle_spindles_input",        # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  5,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["m1ms0","m1ms1","m1ms2","m1ms3","m1ms4"],
            "meshName":                   "muscleSpindleMesh",            # the mesh on which the additional field variables will be defined
            "beforeComputation": [       # transfer/mapping of dofs that will be performed before the computation of the nested solver
              # read spindle stretch (slot m_lda) and communicate to all processes
              {
                "fromConnectorSlot":                "m1lda",
                "toConnectorSlots":                 "m1ms0",
                "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "communicate",        # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "global",
                "toDofNosNumbering":                "global",
                "dofsMapping":
                  {muscle_spindle_dof : [rank_no*variables.n_muscle_spindles + i for rank_no in range(n_ranks)]
                   for i,muscle_spindle_dof in enumerate(muscle_spindle_node_nos)},
                "inputDofs":                        None,
                "outputDofs":                       None,
                "callback":                         None,
                #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              },
              # call callback_muscle_spindles_input
              {
                "fromConnectorSlot":                "m1ms0",
                "toConnectorSlots":                 ["m1ms0","m1ms1","m1ms2","m1ms3","m1ms4"],
                "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        list(range(variables.n_muscle_spindles)),
                "outputDofs":                       [list(range(variables.n_muscle_spindles)) for _ in range(5)],   # [0,1,...,n_muscle_spindles]
                "callback":                         variables.callback_muscle_spindles_input,
                #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
            ],
            "afterComputation":  None,

            # map from motoneuronMesh to stimulated nodes
            "MapDofs": {
              "description":                "motoneurons->stimulated nodes",  # description that will be shown in solver structure visualization
              "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
              "additionalSlotNames":        ["m1mn"],
              "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
              "beforeComputation":          None,

              # mapping from motoneuronMesh which contains on every rank as many nodes as there are motoneurons to the 3D domain
              # map from motoneuronMesh (algebraics) to the fiber meshes (solution)
              "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
                {
                  "fromConnectorSlot":                "m1mn",                # source slot of the dofs mapping
                  "toConnectorSlots":                 "m1vm",                # target slot of the dofs mapping
                  "fromSlotConnectorArrayIndex":      0,
                  "toSlotConnectorArrayIndex":        get_fiber_index_in_motor_unit(fiber_index, motor_unit_no),      # which fiber in this motor unit
                  "mode":                             "localSetIfAboveThreshold",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                  "fromDofNosNumbering":              "local",
                  "toDofNosNumbering":                "global",
                  "dofsMapping":                      {motor_unit_no: stimulation_node_nos},   # map from the motor unit to the stimulated node of the fiber mesh
                  "inputDofs":                        None,                # this option is only needed in mode "callback"
                  "outputDofs":                       None,                # this option is only needed in mode "callback"
                  "callback":                         None,                # this option is only needed in mode "callback"
                  "thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                  "valueToSet":                       variables.vm_value_stimulated,       # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
                }
              for motor_unit_no in range(variables.n_motor_units)
                for fiber_index in range(get_n_fibers_in_motor_unit(motor_unit_no))],     # iterate over all motor units and all fibers in every motor unit


              "Coupling": {
                "timeStepWidth":          variables.dt_bidomain,  # 1e-1
                "logTimeStepWidthAsKey":  "dt_3D",
                "durationLogKey":         "duration_total_muscle1",
                "timeStepOutputInterval": 50,
                "endTime":                1,
                "connectedSlotsTerm1To2": {0:0},  # elasticity: transfer gamma to elasticity, fibers_emg: transfer Vm to StaticBidomainSolver
                "connectedSlotsTerm2To1": [None],   # elasticity: only transfer back geometry (this happens automatically),   fibers_emg: transfer nothing back
                "Term1": {        # monodomain, fibers
                  "MultipleInstances": {
                    "logKey":                     "duration_subdomains_xy_muscle1",
                    "ranksAllComputedInstances":  list(range(n_ranks)),
                    "nInstances":                 variables.n_subdomains_xy,
                    "instances":
                    [{
                      "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
                      "StrangSplitting": {
                        #"numberTimeSteps": 1,
                        "timeStepWidth":          variables.dt_splitting_0D1D,  # 1e-1
                        "logTimeStepWidthAsKey":  "dt_splitting",
                        "durationLogKey":         "duration_monodomain_muscle1",
                        "timeStepOutputInterval": 100,
                        "endTime":                variables.dt_splitting_0D1D,
                        "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion), for elasticity also transfer gamma
                        "connectedSlotsTerm2To1": [0],   # transfer the same back, this avoids data copy

                        "Term1": {      # CellML, i.e. reaction term of Monodomain equation
                          "MultipleInstances": {
                            "logKey":             "duration_subdomains_z_muscle1",
                            "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                            "instances":
                            [{
                              "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                              "Heun" : {
                                "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                                "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                                "durationLogKey":               "duration_0D_muscle1",                           # log key of duration for this solver
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
                                  "setSpecificStatesFunction":              None,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                                  #"setSpecificStatesCallInterval":          2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                                  "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                                  "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                                  "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                                  "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                                  "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                                  "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.

                                  # parameters to the cellml model
                                  "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                                  "mappings":                               variables.muscle1_mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py

                                  "meshName":                               "muscle1_fiber{}".format(fiber_no),                # reference to the fiber mesh
                                  "stimulationLogFilename":                 "out/" + variables.scenario_name + "/stimulation_muscle1.log",                          # a file that will contain the times of stimulations
                                },
                                "OutputWriter" : [
                                  {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/muscle1_0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
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
                                "ImplicitEuler": {
                                  "initialValues":               [],                                      # initial values to be set in the solution vector prior to the first timestep
                                  #"numberTimeSteps":            1,
                                  "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                                  "timeStepWidthRelativeTolerance": 1e-10,                                # tolerance for the time step width, when to rebuild the system matrix
                                  "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                                  "durationLogKey":              "duration_1D_muscle1",                           # log key of duration for this solver
                                  "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep to console
                                  "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                                  "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                                  "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                                  "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                                  "checkForNanInf":              False,                                   # if the solution should be checked for NaN and inf values, this requires a lot of runtimes
                                  "nAdditionalFieldVariables":   2,    # number of additional field variables that should be added and potentially written to output files, these field variables can be used for receiving data from other solvers
                                  "additionalSlotNames":         [],                                      # slot names for the additional field variables
                                  "FiniteElementMethod" : {
                                    "inputMeshIsGlobal":         True,
                                    "meshName":                  "muscle1_fiber{}".format(fiber_no),
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
                              "OutputWriter" : variables.output_writer_fibers_muscle1,
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
                  "OutputSurface": {        # version for fibers_emg_2d_output
                    "OutputWriter": [
                      {"format": "Paraview", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_surface), "filename": "out/" + variables.scenario_name + "/muscle1_surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental",},
                    ] if variables.enable_surface_emg else [],
                    "face":                     ["0+"],              # which faces of the 3D mesh should be written into the 2D mesh
                    "samplingPoints":           None,                # the electrode positions, they are created in the helper.py script
                    "updatePointPositions":     False,               # the electrode points should be initialize in every timestep (set to False for the static case). This makes a difference if the muscle contracts, then True=fixed electrodes, False=electrodes moving with muscle.
                    "filename":                 "out/{}/muscle1_electrodes.csv".format(variables.scenario_name),
                    "enableCsvFile":            False,               # if the values at the sampling points should be written to csv files
                    "enableVtpFile":            False,               # if the values at the sampling points should be written to vtp files
                    "enableGeometryInCsvFile":  False,               # if the csv output file should contain geometry of the electrodes in every time step. This increases the file size and only makes sense if the geometry changed throughout time, i.e. when computing with contraction
                    "enableGeometryFiles":      False,               # if there should be extra files of the locations of the electrodes on every rank
                    "xiTolerance":              0.3,                 # tolerance for element-local coordinates xi, for finding electrode positions inside the elements. Increase or decrease this numbers if not all electrode points are found.

                    "StaticBidomainSolver": {
                      "timeStepWidth":          variables.dt_bidomain,
                      "timeStepOutputInterval": 50,
                      "durationLogKey":         "duration_bidomain_muscle1",
                      "solverName":             "muscularEMGSolver",
                      "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
                      "enableJacobianConditionNumber": False,        # if set to true, estimate the condition number of the jacobian of the element-coordinate-to-world-frame mapping in every element and output it in the output writer
                      "slotNames":              [],
                      "nAdditionalFieldVariables": 0,

                      "PotentialFlow": {
                        "FiniteElementMethod" : {
                          "meshName":           "muscle1Mesh",
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
                          "meshName":           "muscle1Mesh",
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
                      "OutputWriter" : variables.output_writer_emg_muscle1,
                    }
                  }
                }
              }
            }
          }
        },
        "Term5": {
          # muscle2: bidoamin + 1D monodomain + 0D
          "Coupling": {
            "timeStepWidth":          variables.dt_bidomain,  # 1e-1
            "logTimeStepWidthAsKey":  "dt_3D",
            "durationLogKey":         "duration_total_muscle2",
            "timeStepOutputInterval": 50,
            "endTime":                1,
            "connectedSlotsTerm1To2": {0:0},  # elasticity: transfer gamma to elasticity, fibers_emg: transfer Vm to StaticBidomainSolver
            "connectedSlotsTerm2To1": [None],   # elasticity: only transfer back geometry (this happens automatically),   fibers_emg: transfer nothing back
            "Term1": {        # monodomain, fibers
              "MultipleInstances": {
                "logKey":                     "duration_subdomains_xy_muscle2",
                "ranksAllComputedInstances":  list(range(n_ranks)),
                "nInstances":                 variables.n_subdomains_xy,
                "instances": 
                [{
                  "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
                  "StrangSplitting": {
                    #"numberTimeSteps": 1,
                    "timeStepWidth":          variables.dt_splitting_0D1D,  # 1e-1
                    "logTimeStepWidthAsKey":  "dt_splitting",
                    "durationLogKey":         "duration_monodomain_muscle2",
                    "timeStepOutputInterval": 100,
                    "endTime":                variables.dt_splitting_0D1D,
                    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion), for elasticity also transfer gamma
                    "connectedSlotsTerm2To1": [0],   # transfer the same back, this avoids data copy

                    "Term1": {      # CellML, i.e. reaction term of Monodomain equation
                      "MultipleInstances": {
                        "logKey":             "duration_subdomains_z_muscle2",
                        "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                        "instances": 
                        [{
                          "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                          "Heun" : {
                        "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                        "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                        "durationLogKey":               "duration_0D_muscle2",                           # log key of duration for this solver
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
                          "setSpecificStatesFunction":              None,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                          #"setSpecificStatesCallInterval":          2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                          "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                          "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                          "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                          "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                          "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                          "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                          
                          # parameters to the cellml model
                          "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                          "mappings":                               variables.muscle2_mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                          
                          "meshName":                               "muscle2_fiber{}".format(fiber_no),                # reference to the fiber mesh
                          "stimulationLogFilename":                 "out/" + variables.scenario_name + "/stimulation_muscle2.log",                          # a file that will contain the times of stimulations
                        },      
                        "OutputWriter" : [
                          {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/muscle2_0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
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
                          "ImplicitEuler": {
                        "initialValues":               [],                                      # initial values to be set in the solution vector prior to the first timestep
                        #"numberTimeSteps":            1,
                        "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                        "timeStepWidthRelativeTolerance": 1e-10,                                # tolerance for the time step width, when to rebuild the system matrix
                        "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                        "durationLogKey":              "duration_1D_muscle2",                           # log key of duration for this solver
                        "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep to console
                        "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                        "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                        "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                        "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                        "checkForNanInf":              False,                                   # if the solution should be checked for NaN and inf values, this requires a lot of runtimes
                        "nAdditionalFieldVariables":   2,    # number of additional field variables that should be added and potentially written to output files, these field variables can be used for receiving data from other solvers
                        "additionalSlotNames":         [],                                      # slot names for the additional field variables
                        "FiniteElementMethod" : { 
                          "inputMeshIsGlobal":         True,
                          "meshName":                  "muscle2_fiber{}".format(fiber_no),
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
                        "OutputWriter" : variables.output_writer_fibers_muscle2,
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
              "OutputSurface": {        # version for fibers_emg_2d_output
                "OutputWriter": [
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_bidomain*variables.output_timestep_surface), "filename": "out/" + variables.scenario_name + "/muscle2_surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental",},
                ] if variables.enable_surface_emg else [],
                "face":                     ["0+"],              # which faces of the 3D mesh should be written into the 2D mesh
                "samplingPoints":           None,                # the electrode positions, they are created in the helper.py script
                "updatePointPositions":     False,               # the electrode points should be initialize in every timestep (set to False for the static case). This makes a difference if the muscle contracts, then True=fixed electrodes, False=electrodes moving with muscle.
                "filename":                 "out/{}/muscle2_electrodes.csv".format(variables.scenario_name),
                "enableCsvFile":            False,               # if the values at the sampling points should be written to csv files
                "enableVtpFile":            False,               # if the values at the sampling points should be written to vtp files
                "enableGeometryInCsvFile":  False,               # if the csv output file should contain geometry of the electrodes in every time step. This increases the file size and only makes sense if the geometry changed throughout time, i.e. when computing with contraction
                "enableGeometryFiles":      False,               # if there should be extra files of the locations of the electrodes on every rank
                "xiTolerance":              0.3,                 # tolerance for element-local coordinates xi, for finding electrode positions inside the elements. Increase or decrease this numbers if not all electrode points are found.
                
                "StaticBidomainSolver": {
                  "timeStepWidth":          variables.dt_bidomain,
                  "timeStepOutputInterval": 50,
                  "durationLogKey":         "duration_bidomain_muscle2",
                  "solverName":             "muscularEMGSolver",
                  "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
                  "enableJacobianConditionNumber": False,        # if set to true, estimate the condition number of the jacobian of the element-coordinate-to-world-frame mapping in every element and output it in the output writer
                  "slotNames":              [],
                  "nAdditionalFieldVariables": 0,

                  "PotentialFlow": {
                    "FiniteElementMethod" : {
                      "meshName":           "muscle2Mesh",
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
                      "meshName":           "muscle2Mesh",
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
                  "OutputWriter" : variables.output_writer_emg_muscle2,
                }
              }
            }
          }
        }
      }
    },
    "Term2": {
      # 2x mechanics solver: one for each muscle
      "Coupling": {
        "timeStepWidth":          variables.dt_elasticity,
        "logTimeStepWidthAsKey":  "dt_elasticity",
        "durationLogKey":         "duration_elasticity",
        "timeStepOutputInterval": 1,
        "endTime":                variables.end_time,
        "connectedSlotsTerm1To2": None,       # connect lambda to slot 0 and gamma to slot 2
        "connectedSlotsTerm2To1": None,       # transfer nothing back

        "Term1": {
          "MuscleContractionSolver": {
            "numberTimeSteps":              1,                         # only use 1 timestep per interval
            "timeStepOutputInterval":       1,
            "Pmax":                         variables.Pmax,            # maximum PK2 active stress
            "enableForceLengthRelation":    True,                      # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
            "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
            "slotNames":                    ["m1lda", "m1ldot", "m1g_in", "m1T", "m1ux", "m1uy", "m1uz"],  # slot names of the data connector slots: lambda, lambdaDot, gamma, traction
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + variables.scenario_name + "/muscle1_contraction", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
            ],
            "mapGeometryToMeshes":          ["muscle1Mesh"] + [key for key in fiber_meshes.keys() if "muscle1_fiber" in key],    # the mesh names of the meshes that will get the geometry transferred
            "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
            "dynamic":                      variables.dynamic,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
            
            # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
            "DynamicHyperelasticitySolver": {
              "timeStepWidth":              variables.dt_elasticity,           # time step width 
              "durationLogKey":             "muscle1_duration_mechanics",               # key to find duration of this solver in the log file
              "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
              
              "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
              "density":                    variables.rho,             # density of the material
              "dampingFactor":              variables.damping_factor,  # factor for velocity dependent damping
              "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
              "residualNormLogFilename":    "out/"+variables.scenario_name+"/muscle1_log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
              "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
              "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
                
              "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
              # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
              
              # mesh
              "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
              "meshName":                   "muscle1Mesh_quadratic",       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
              "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction
              "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
        
              # solving
              "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
              #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
              "loadFactors":                [],                        # no load factors, solve problem directly
              "loadFactorGiveUpThreshold":  0.25,                       # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
              "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
              "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
              
              # boundary and initial conditions
              "dirichletBoundaryConditions": variables.muscle1_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
              "neumannBoundaryConditions":   variables.muscle1_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
              "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
              "updateDirichletBoundaryConditionsFunction": None,                  # muscle1_update_dirichlet_boundary_conditions_helper, function that updates the dirichlet BCs while the simulation is running
              "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
              "updateNeumannBoundaryConditionsFunction":   muscle1_update_neumann_boundary_conditions_helper,                    # function that updates the Neumann BCs while the simulation is running
              "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step

              
              "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
              "constantBodyForce":           variables.main_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
              
              "dirichletOutputFilename":     "out/"+variables.scenario_name+"/muscle1_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
              "totalForceLogFilename":       "out/"+variables.scenario_name+"/muscle1_tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
              "totalForceLogOutputInterval":       10,                                  # output interval when to write the totalForceLog file

              # define which file formats should be written
              # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
              "OutputWriter" : [
                
                # Paraview files
                {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_displacements", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                
                # Python callback function "postprocess"
                # responsible to model the tendon
                {"format": "PythonCallback", "outputInterval": 1, "callback": variables.muscle1_postprocess, "onlyNodalValues":True, "filename": "", "fileNumbering":'incremental'},
              ],
              # 2. additional output writer that writes also the hydrostatic pressure
              "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
                "OutputWriter" : [
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_pressure", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ]
              },
              # 3. additional output writer that writes virtual work terms
              "dynamic": {    # output of the dynamic solver, has additional virtual work values 
                "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                  {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/muscle1_dynamic", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ],
              },
              # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
              "LoadIncrements": {   
                "OutputWriter" : [
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle1_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ]
              }
            }
          }
        },
        "Term2": {
          "MuscleContractionSolver": {
            "numberTimeSteps":              1,                         # only use 1 timestep per interval
            "timeStepOutputInterval":       1,
            "Pmax":                         variables.Pmax,            # maximum PK2 active stress
            "enableForceLengthRelation":    True,                      # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
            "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
            "slotNames":                    ["m2lda", "m2ldot", "m2g_in", "m2T", "m2ux", "m2uy", "m2uz"],  # slot names of the data connector slots: lambda, lambdaDot, gamma, traction
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + variables.scenario_name + "/muscle2_contraction", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
            ],
            "mapGeometryToMeshes":          ["muscle2Mesh"] + [key for key in fiber_meshes.keys() if "muscle2_fiber" in key],    # the mesh names of the meshes that will get the geometry transferred
            "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
            "dynamic":                      variables.dynamic,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
            
            # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
            "DynamicHyperelasticitySolver": {
              "timeStepWidth":              variables.dt_elasticity,           # time step width 
              "durationLogKey":             "muscle2_duration_mechanics",               # key to find duration of this solver in the log file
              "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
              
              "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
              "density":                    variables.rho,             # density of the material
              "dampingFactor":              variables.damping_factor,  # factor for velocity dependent damping
              "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
              "residualNormLogFilename":    "out/"+variables.scenario_name+"/muscle2_log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
              "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
              "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
                
              "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
              # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
              
              # mesh
              "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
              "meshName":                   "muscle2Mesh_quadratic",       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
              "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction
              "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
        
              # solving
              "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
              #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
              "loadFactors":                [],                        # no load factors, solve problem directly
              "loadFactorGiveUpThreshold":  0.25,                       # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
              "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
              "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
              
              # boundary and initial conditions
              "dirichletBoundaryConditions": variables.muscle2_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
              "neumannBoundaryConditions":   variables.muscle2_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
              "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
              "updateDirichletBoundaryConditionsFunction": None,                  # muscle2_update_dirichlet_boundary_conditions_helper, function that updates the dirichlet BCs while the simulation is running
              "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
              "updateNeumannBoundaryConditionsFunction":   muscle2_update_neumann_boundary_conditions_helper,                    # function that updates the Neumann BCs while the simulation is running
              "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step

              
              "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
              "constantBodyForce":           variables.main_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
              
              "dirichletOutputFilename":     "out/"+variables.scenario_name+"/muscle2_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
              "totalForceLogFilename":       "out/"+variables.scenario_name+"/muscle2_tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
              "totalForceLogOutputInterval":       10,                                  # output interval when to write the totalForceLog file

              # define which file formats should be written
              # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
              "OutputWriter" : [
                
                # Paraview files
                {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle2_displacements", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                
                # Python callback function "postprocess"
                # responsible to model the tendon
                {"format": "PythonCallback", "outputInterval": 1, "callback": variables.muscle2_postprocess, "onlyNodalValues":True, "filename": "", "fileNumbering":'incremental'},
              ],
              # 2. additional output writer that writes also the hydrostatic pressure
              "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
                "OutputWriter" : [
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle2_pressure", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ]
              },
              # 3. additional output writer that writes virtual work terms
              "dynamic": {    # output of the dynamic solver, has additional virtual work values 
                "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                  {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/muscle2_dynamic", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle2_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ],
              },
              # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
              "LoadIncrements": {   
                "OutputWriter" : [
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/muscle2_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ]
              }
            }
          }
        }
      }
    }
  }
}

# stop   timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))

