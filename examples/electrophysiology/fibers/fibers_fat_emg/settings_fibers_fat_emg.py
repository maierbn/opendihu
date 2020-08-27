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
  variables_file = sys.argv[0]
  variables_module = variables_file[0:variables_file.find(".py")]
  
  if rank_no == 0:
    print("Loading variables from {}.".format(variables_file))
    
  custom_variables = importlib.import_module(variables_module)
  variables.__dict__.update(custom_variables.__dict__)
  sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed
else:
  if rank_no == 0:
    print("Error: There is no variables file, e.g:\n ./fibers_fat_emg ../settings_fibers_fat_emg.py ramp.py\n")
  exit(0)

# -------------- begin user parameters ----------------

# -------------- end user parameters ----------------

# define command line arguments
parser = argparse.ArgumentParser(description='fibers_fat_emg')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--diffusion_solver_type',               help='The solver for the diffusion.',               default=variables.diffusion_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--diffusion_preconditioner_type',       help='The preconditioner for the diffusion.',       default=variables.diffusion_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).', default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',  default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_type',                     help='The solver for the static bidomain.',         default=variables.emg_solver_type)
#parser.add_argument('--emg_solver_type',                    help='The solver for the static bidomain.',         default=variables.emg_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--emg_preconditioner_type',             help='The preconditioner for the static bidomain.', default=variables.emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_initial_guess_nonzero',           help='If the initial guess for the emg linear system should be set to the previous solution.', default=variables.emg_initial_guess_nonzero, action='store_true')
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',           type=float, default=variables.output_timestep)
parser.add_argument('--dt_0D',                               help='The timestep for the 0D model.',              type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D',                               help='The timestep for the 1D model.',              type=float, default=variables.dt_1D)
parser.add_argument('--dt_splitting',                        help='The timestep for the splitting.',             type=float, default=variables.dt_splitting)
parser.add_argument('--dt_3D',                               help='The timestep for the 3D model, either bidomain or mechanics.', type=float, default=variables.dt_3D)
parser.add_argument('--disable_firing_output',               help='Disables the initial list of fiber firings.', default=variables.disable_firing_output, action='store_true')
parser.add_argument('--fast_monodomain_solver_optimizations',help='Enable the optimizations for fibers.',        default=variables.fast_monodomain_solver_optimizations, action='store_true')
parser.add_argument('--use_vc',                              help='If the Vc optimization type should be used for cellml adapter.',  default=variables.use_vc, action='store_true')
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
    for j in [1]:     # partitioning in y direction is always set to 1 because the composite mesh of muscle volume and surface needs to have all processes use all submeshes and the fat layer mesh would not have subdomains for all processes if both the partitioning in x and y direction would be !=1
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
  print("dt_0D:           {:0.0e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.0e}, potential_flow_solver_type: {}".format(variables.dt_1D, variables.potential_flow_solver_type))
  print("dt_splitting:    {:0.0e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.dt_splitting, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("fast_monodomain_solver_optimizations: {}, use_vc: {}".format(variables.fast_monodomain_solver_optimizations, variables.use_vc))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
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

# Function to postprocess the output
# This function gets periodically called by the running simulation. 
# It provides all current variables for each node
def postprocess(result):
  result = result[0]
  current_time = result["currentTime"]
  field_variables = result["data"]
  
  # result contains the most recent data
  # field_variables[0] is geometry (3 components "x","y","z")
  # field_variables[1] is fiberDirection (3 components)
  # field_variables[2] is phi_e (1 component)
  # field_variables[3] is Vm (1 component)
  # field_variables[4] is transmembraneFlow (1 component)
  # field_variables[5] is flowPotential (1 component)
  
  # get all emg values
  phi_e_values = field_variables[2]["components"][0]["values"]
  geometry_x_values = field_variables[0]["components"][0]["values"]
  geometry_y_values = field_variables[0]["components"][1]["values"]
  geometry_z_values = field_variables[0]["components"][2]["values"]
  
  # select nodes of the fat layer mesh, z = direction along muscle, x = across muscle (i.e. in x and y direction)
  x_begin = 2           # first index to select
  x_end = x_begin + 13   # one after last index to select
  x_step = 2            # stride which node to select
  
  # in fiber direction
  z_begin = 4
  z_end = z_begin + 20
  z_step = 2
  
  # get helper variables, dimensions of the fat layer mesh
  n_points_local_x = variables.fat_mesh_n_points_local[0]
  n_points_local_y = variables.fat_mesh_n_points_local[1]
  n_points_local_z = variables.fat_mesh_n_points_local[2]
  
  n_points_global_x = variables.fat_mesh_n_points_global[0]
  n_points_global_y = variables.fat_mesh_n_points_global[1]
  n_points_global_z = variables.fat_mesh_n_points_global[2]
  
  index_offset_x = variables.fat_mesh_index_offset[0]
  index_offset_y = variables.fat_mesh_index_offset[1]
  index_offset_z = variables.fat_mesh_index_offset[2]
  
  #print("current_time: {}, index offset: {},{},{}".format(current_time, index_offset_x, index_offset_y, index_offset_z))
  
  # loop over the global indices and only continue if they are on the local domain
  for z_index in range(z_begin, z_end, z_step):
    if not (index_offset_z <= z_index < index_offset_z+n_points_local_z):
      continue
    for x_index in range(x_begin, x_end, x_step):
      if not (index_offset_x <= x_index < index_offset_x+n_points_local_x):
        continue
      
      # here, (x_index,z_index) are the global indices of a point that is in the local domain
      
      # get the local coordinates for (x_index,z_index)
      k_local = z_index - index_offset_z
      j_local = n_points_local_y-1        # select top layer of fat mesh
      i_local = x_index - index_offset_x
      
      index = k_local*n_points_local_y*n_points_local_x + j_local*n_points_local_x + i_local
      
      # output message if index is out of bounds
      if index >= len(phi_e_values):
        print("{}: local index {} ({},{},{}) is >= size {} ({} x {} x {}) ".format(rank_no, index, i_local, j_local, k_local, len(phi_e_values), n_points_local_x, n_points_local_y, n_points_local_z))
      
      # get the value
      phi_e_value = phi_e_values[index]
      
      # save to file
      filename = "out/emg_{:02}_{:02}.csv".format(x_index,z_index)
      
      # clear file at the beginning of the simulation
      if current_time <= variables.output_timestep_electrodes + 1e-5:
        with open(filename,"w") as f:
          f.write("time;phi_e;\n".format(current_time,phi_e_value))
          
        filename_points = "out/points.{}.csv".format(rank_no)
        with open(filename_points,"a") as f:
          f.write("{};{};{}\n".format(geometry_x_values[index], geometry_y_values[index], geometry_z_values[index]))
      
      # append line with current time and EMG value
      with open(filename,"a") as f:
        f.write("{};{};\n".format(current_time,phi_e_value))

# define the config dict
config = {
  "scenarioName":                  variables.scenario_name,
  "mappingsBetweenMeshesLogFile":  "out/mappings_between_meshes.txt",
  "logFormat":                     "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":    "out/solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "meta": {                 # additional fields that will appear in the log
    "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": variables.mappings_between_meshes,
  "Solvers": {
    "diffusionTermSolver": {# solver for the implicit timestepping scheme of the diffusion time step
      "maxIterations":      1e4,
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    },
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
  "Coupling": {
    "timeStepWidth":          variables.dt_3D,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_3D",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 10,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": [0],          # transfer Vm to StaticBidomainSolver
    "connectedSlotsTerm2To1": [None],       # transfer nothing back
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
            "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
            "connectedSlotsTerm2To1": [0],   # transfer the same back, this avoids data copy

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
                    "checkForNanInf":               True,                                    # abort execution if the solution contains nan or inf values
                    "nAdditionalFieldVariables":    0,                                       # number of additional field variables
                    "additionalSlotNames":          [],                                      # names for the additional slots
                      
                    "CellML" : {
                      "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                      #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                      "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                      "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                      
                      # optimization parameters
                      "optimizationType":                       "vc" if variables.use_vc else "simd",           # "vc", "simd", "openmp" type of generated optimizated source file
                      "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                      "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                      "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                      
                      # stimulation callbacks
                      #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                      #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                      #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                      #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                      "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                      "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                      "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                      "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                      "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                      
                      # parameters to the cellml model
                      "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                      "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                      
                      "meshName":                               "MeshFiber_{}".format(fiber_no),                # reference to the fiber mesh
                      "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
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
                  "ranks":                         list(range(variables.n_subdomains_z)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "ImplicitEuler" : {
                    "initialValues":               [],                                      # no initial values are given
                    #"numberTimeSteps":            1,
                    "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                    "timeStepWidthRelativeTolerance": 1e-10,
                    "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                    "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                    "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep
                    "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                    "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                    "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                    "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                    "nAdditionalFieldVariables":   0,
                    "additionalSlotNames":         [],
                    "checkForNanInf":              False,
                    
                    "FiniteElementMethod" : {
                      "inputMeshIsGlobal":         True,
                      "meshName":                  "MeshFiber_{}".format(fiber_no),
                      "solverName":                "diffusionTermSolver",
                      "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                      "slotName":                  "vm",
                    },
                    "OutputWriter" : [
                      #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                      #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                      #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                      #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
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
    },
    "Term2": {        # Bidomain, EMG
      "OutputSurface": {
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_surface), "filename": "out/" + variables.scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
        ],
        #"face":                    ["1+","0+"],         # which faces of the 3D mesh should be written into the 2D mesh
        "face":                     ["1+"],              # which faces of the 3D mesh should be written into the 2D mesh
        "samplingPoints":           variables.hdemg_electrode_positions,    # the electrode positions, they are created in the helper.py script
        "updatePointPositions":     False,               # the electrode points should be initialize in every timestep (set to False for the static case). This makes a difference if the muscle contracts, then True=fixed electrodes, False=electrodes moving with muscle.
        "filename":                 "out/{}/electrodes.csv".format(variables.scenario_name),
        "xiTolerance":              0.3,                 # tolerance for element-local coordinates xi, for finding electrode positions inside the elements. Increase or decrease this numbers if not all electrode points are found.
        "StaticBidomainSolver": {             # solves Bidomain equation: K(sigma_i) Vm + K(sigma_i+sigma_e) phi_e = 0   => K(sigma_i+sigma_e) phi_e = -K(sigma_i) Vm
          "numberTimeSteps":        1,
          "timeStepOutputInterval": 50,
          "durationLogKey":         "duration_bidomain",
          "solverName":             "muscularEMGSolver",
          "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
          "slotNames":              [],
          "PotentialFlow": {
            "FiniteElementMethod" : {  
              "meshName":           ["3Dmesh","3DFatMesh"],
              "solverName":         "potentialFlowSolver",
              "prefactor":          1.0,
              "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
              "dirichletOutputFilename":     None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "neumannBoundaryConditions":   [],
              "inputMeshIsGlobal":  True,
              "slotName":           "",
            },
          },
          "Activation": {
            "FiniteElementMethod" : {  
              "meshName":           ["3Dmesh","3DFatMesh"],       # composite mesh that consists of the muscle value mesh and the fat layer mesh
              "solverName":         "muscularEMGSolver",
              "prefactor":          1.0,
              "inputMeshIsGlobal":  True,
              "dirichletBoundaryConditions": {},
              "dirichletOutputFilename":     None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
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
          "OutputWriter" : variables.output_writer_emg + [
            {"format": "PythonCallback", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_electrodes), "onlyNodalValues":True, "filename": "", "callback": postprocess, "fileNumbering": "incremental"}
          ],
        }
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
