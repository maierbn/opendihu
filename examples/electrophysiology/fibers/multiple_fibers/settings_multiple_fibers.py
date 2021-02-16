# Multiple fibers from biceps geometry.
# This example has MegaMol integration but also outputs Paraview files.
# The number of fibers depends on the number of processes.
#
# arguments:  [<n_processes_per_fiber> [<scenario_name>]]
#
# E.g. to have 2 fibers with 2 processes, each:
# mpirun -n 4 ./multiple_fibers ../settings_multiple_fibers.py 2
#
# It is not possible with this example to have cube-shaped partitions because of the solver structure. 
# In order to have a process compute multiple fibers but only a part of them, use the multiple_fibers_cubes_partitioning example.
# Compare  multiple_fibers_cubes_partitioning/src/multiple_fibers.cpp  with multiple_fibers/src/multiple_fibers.cpp to get the difference.

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import sys
import argparse

import py_reader

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm

megamol_output_timestep = 5e-2    # timestep frequency for megamol

# default parameter values
class Variables:
    scenario_name = ''
    n_subdomains_z = 1
    diffusion_solver_type = 'cg'
    diffusion_preconditioner_type = 'none'
    diffusion_solver_maxit = 1e4
    diffusion_solver_reltol = 1e-10

    #fiber_file = "../../../input/laplace3d_structured_quadratic"
    fiber_file = "../../../input/laplace3d_structured_linear"
    #fiber_file = "../../../input1000/laplace3d_structured_quadratic"

    # input files
    #cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
    #cellml_file = "../../../input/shorten.cpp"
    cellml_file = "../../../input/hodgkin_huxley_1952.c"

    fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
    #firing_times_file = "../../../input/MU_firing_times_real.txt"
    firing_times_file = "../../../input/MU_firing_times_immediately.txt"
    stimulation_frequency = 9.0       # stimulations per ms => each stimulus takes 1ms/stimulation_frequency

    end_time = 20.0
    output_timestep = 4e-1            # timestep for output files
    output_interval_0D = 0            # output iterval vor 0D model (0 = 'use output timestep')
    output_interval_1D = 0            # output iterval vor 1D model (0 = 'use output timestep')
    dt_1D = 1e-3                      # timestep width of diffusion
    dt_0D = 3e-3                      # timestep width of ODEs
    dt_splitting = 3e-3               # overall timestep width of splitting
    nt_1D = 0                         # number of timesteps per 0D splitting substep (0 = 'use time step width')
    nt_0D = 0                         # number of timesteps per 1D splitting substep (0 = 'use time step width')
    nt_splitting = 0                  # number of splitting timesteps (0 = 'use time step width')

    time_step_width_rel_tol_1D = 1e-10

    initial_value_file = None         # python file which contians the initial values. VTK only stores 32bit floats.
    disable_firing = np.infty         # time after which we disable the firing

    outfile_0D = "out/fibre_{i}_0D"
    outfile_1D = "out/fibre_{i}_1D"

variables = Variables()

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

parser = argparse.ArgumentParser(description='fibers_emg, All times are in 1ms.')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--diffusion_solver_type',               help='The solver for the diffusion.',               default=variables.diffusion_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--diffusion_preconditioner_type',       help='The preconditioner for the diffusion.',       default=variables.diffusion_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--diffusion_solver_reltol',             help='Ralative tolerance for diffusion solver',     type=float, default=variables.diffusion_solver_reltol)
parser.add_argument('--diffusion_solver_maxit',              help='Maximum number of iterations for diffusion solver', type=int, default=variables.diffusion_solver_maxit)
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--cellml_file',                         help='The filename of the file that contains the cellml model.', default=variables.cellml_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--stimulation_frequency',               help='Stimulations per ms. Each stimulation corresponds to one line in the firing_times_file.', default=variables.stimulation_frequency)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',           type=float, default=variables.output_timestep)
parser.add_argument('--output_interval_0D',                  help='The interval for writing 0D outputs. Overrides `--output_timestep` for 0D model.', type=int, default=variables.output_interval_0D)
parser.add_argument('--output_interval_1D',                  help='The interval for writing 1D outputs. Overrides `--output_timestep` for 1D model.', type=int, default=variables.output_interval_1D)
parser.add_argument('--dt_0D',                               help='The timestep for the 0D model.',              type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D',                               help='The timestep for the 1D model.',              type=float, default=variables.dt_1D)
parser.add_argument('--dt_splitting',                        help='The timestep for the splitting.',             type=float, default=variables.dt_splitting)
parser.add_argument('--nt_0D',                               help='The number of timesteps for the 0D model. Overrides `--dt_0D`.', type=int, default=variables.nt_0D)
parser.add_argument('--nt_1D',                               help='The number of timesteps for the 1D model. Overrides `--dt_1D`.', type=int, default=variables.nt_1D)
parser.add_argument('--nt_splitting',                        help='The number of splitting timesteps to reach `--tend`. Overrides `--dt_splitting`.', type=int, default=variables.nt_splitting)
parser.add_argument('--time_step_width_rel_tol_1D',          help='The relative tolerance for the time step size used to recompute the system matrix.', type=float, default=variables.time_step_width_rel_tol_1D)
parser.add_argument('--initial_value_file',                  help='Initial value for V,m,h,n. Only python files are pupported.', default=variables.initial_value_file)
parser.add_argument('--disable_firing',                      help='Disable stimulus after certain time. Useful in combination with --initial_values', type=float, default=variables.disable_firing)
parser.add_argument('--outfile_0D',                          help='Output file name for 0D time steps. Use {i} for fiber index. Set to empty to disable output', default=variables.outfile_0D)
parser.add_argument('--outfile_1D',                          help='Output file name for 1D time steps. Use {i} for fiber index. Set to empty to disable output', default=variables.outfile_1D)
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

if variables.n_subdomains is not None:
    variables.n_subdomains_x = variables.n_subdomains[0]
    variables.n_subdomains_y = variables.n_subdomains[1]
    variables.n_subdomains_z = variables.n_subdomains[2]
n_processes_per_fiber = variables.n_subdomains_z

if rank_no == 0:
  print("n_processes_per_fiber: {}".format(n_processes_per_fiber))
  print("scenario_name: {}".format(variables.scenario_name))

#print("rank: {}/{}".format(rank_no,n_ranks))

# set values for cellml model
if "shorten" in variables.cellml_file:
  parameters_used_as_algebraic = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 400.

elif "hodgkin_huxley" in variables.cellml_file:
  parameters_used_as_algebraic = []
  parameters_used_as_constant = [2]
  parameters_initial_values = [0.0]
  nodal_stimulation_current = 40.

def get_motor_unit_no(fibre_no):
  return int(fibre_distribution[fibre_no % len(fibre_distribution)]-1)

def fiber_gets_stimulated(fibre_no, frequency, current_time):
  # frequency is number of (potential) stimulations per ms / how many rows are to be read
  # if the fiber is stumulated depends on the value in the "firing times tile"

  # determine motor unit
  mu_no = (int)(get_motor_unit_no(fibre_no)*0.8)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  return firing_times[index % firing_times.shape[0], mu_no] == 1 # each column correspond to one motor unit

def set_parameters_null(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  pass
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fibre_no, variables.stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)
  
  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = nodal_stimulation_current
  else:
    stimulation_current = 0.
  
  first_dof_global = dof_nos_global[0]
  last_dof_global = dof_nos_global[-1]
    
  for node_no_global in nodes_to_stimulate_global:
    if first_dof_global <= node_no_global <= last_dof_global:
      # get local no for global no (1D)
      dof_no_local = node_no_global - first_dof_global
      parameters[dof_no_local] = stimulation_current
 
      #print("       {}: set stimulation for local dof {}".format(rank_no, dof_no_local))
  
  #print("       {}: setParameters at timestep {}, t={}, n_nodes_global={}, range: [{},{}], fibre no {}, MU {}, stimulated: {}".\
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fibre_no, get_motor_unit_no(fibre_no), is_fiber_gets_stimulated))
    
  #wait = input("Press any key to continue...")
    
# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fibre_no, variables.stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(10):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
  
  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = 40.
  else:
    stimulation_current = 0.

  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index)


if variables.initial_value_file is not None:
    def extract4(file):
        if file.endswith('.py'):
            data = py_reader.load_data([file])
            data = data[0]['data']
            solution  = next(filter(lambda d: d['name'] == 'solution', data))
            componentX = lambda x: next(filter(lambda d: d['name'] == str(x), solution['components']))
            channel_names = ['membrane/V', 'sodium_channel_m_gate/m', 'sodium_channel_h_gate/h', 'potassium_channel_n_gate/n']
            return np.vstack([componentX(i)['values'] for i in channel_names]).T
            raise "FileType not understood: "+file
    print("load initial state from "+variables.initial_value_file+" for V,m,h,n")
    data = extract4(variables.initial_value_file)
    initial_values = [list(states) for states in data]; # [(V m h n)...]
    initial_values_cellml = "undefined"
else:
    # use the initial values defined in the CellML model
    initial_values = []
    initial_values_cellml = "CellML"

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fibre_no):
  if current_time < variables.disable_firing:
      # determine if fibre gets stimulated at the current time
      is_fiber_gets_stimulated = fiber_gets_stimulated(fibre_no, variables.stimulation_frequency, current_time)

      if is_fiber_gets_stimulated:
        # determine nodes to stimulate (center node, left and right neighbour)
        innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
        innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
        nodes_to_stimulate_global = [innervation_node_global]

        for node_no_global in nodes_to_stimulate_global:
          states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)

def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  pass
    
def get_instance_config(i):

  # set ranks list containing the rank nos for fiber i  
  ranks = []
  for j in range(n_processes_per_fiber):
    ranks.append(n_processes_per_fiber*i + j)

  if variables.output_interval_0D == 0:
    output_interval_0D = int(1./variables.dt_0D*variables.output_timestep)
  else:
    output_interval_0D = variables.output_interval_0D
  if variables.output_interval_1D == 0:
    output_interval_1D = int(1./variables.dt_1D*variables.output_timestep)
  else:
    output_interval_1D = variables.output_interval_1D
  bc = {0: -75, -1: -75}
  instance_config = {
    "ranks": ranks,
    "StrangSplitting": {
      "timeStepWidth":          variables.dt_splitting,  # 1e-1
      "numberTimeSteps":        variables.nt_splitting,
      "logTimeStepWidthAsKey":  "dt_splitting",
      "durationLogKey":         "duration_total",
      "timeStepOutputInterval": 1000,
      "endTime":                variables.end_time,
      "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
      "connectedSlotsTerm2To1": [0],   # transfer the same back

      "Term1": {      # CellML
        "Heun" : {
          "timeStepWidth":               variables.dt_0D,  # 5e-5
          "numberTimeSteps":             variables.nt_0D,
          "logTimeStepWidthAsKey":       "dt_0D",
          "durationLogKey":              "duration_0D",
          "initialValues":               initial_values,
          "timeStepOutputInterval":      1e4,
          "inputMeshIsGlobal":           True,
          "dirichletBoundaryConditions": {},
          "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
          "nAdditionalFieldVariables":   0,
          "additionalSlotNames":         [],
          "checkForNanInf":              False,
            
          "CellML" : {
            "statesInitialValues":                    initial_values_cellml,
            "modelFilename": variables.cellml_file,   # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            "compilerFlags":                          "-fPIC -O3 -march=native -shared ",
            "optimizationType":                       "vc",     # "vc", "simd", "openmp" type of generated optimizated source file
            "initializeStatesToEquilibrium":          False,
            "approximateExponentialFunction":         False,    # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
            "maximumNumberOfThreads":                 0,        # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
            #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
            #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
            "setSpecificStatesCallInterval": int(1./variables.stimulation_frequency/variables.dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setSpecificStatesCallFrequency":         0,        # set_specific_states should be called variables.stimulation_frequency times per ms, 0 means disabled
            "setSpecificStatesFrequencyJitter":       0,        # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
            "setSpecificStatesRepeatAfterFirstCall":  0,        # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
            "setSpecificStatesCallEnableBegin":       0,        # [ms] first time when to call setSpecificStates
            
            "additionalArgument": i,
            
            "algebraicsForTransfer":  [],            # which algebraic values to use in further computation
            "parametersForTransfer":  [],
            "statesForTransfer":   0,                   # which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
                      
            "parametersUsedAsAlgebraic": parameters_used_as_algebraic,  #[32],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFiber"+str(i),
            "stimulationLogFilename": "out/stimulation.log",
            
          },
          "OutputWriter" : [
            {"format": "Paraview",   "outputInterval": output_interval_0D, "filename": variables.outfile_0D.format(i=i), "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering":"incremental"},
            {"format": "PythonFile", "outputInterval": output_interval_0D, "filename": variables.outfile_0D.format(i=i), "binary": True, "onlyNodalValues":False, "fileNumbering":"incremental"}, # also derivatives for hermite
            # {"format":"PythonCallback", "outputInterval": output_interval_0D, "callback": lambda x: print("writing fiber {ix:5}: 0D for time step {timeStepNo:5} at time {currentTime}".format(ix=i, **x[i]))},
          ] if variables.outfile_0D != '' else []
        },
      },
      "Term2": {     # Diffusion
        "ImplicitEuler" : {
          "initialValues":               [],
          "timeStepWidth":               variables.dt_1D,  # 1e-5
          "numberTimeSteps":             variables.nt_1D,
          "timeStepWidthRelativeTolerance": variables.time_step_width_rel_tol_1D,
          "timeStepWidthRelativeToleranceAsKey": "timeStepRelTol_1D",
          "durationInitTimeStepLogKey":  "duration_init_1D",
          "logTimeStepWidthAsKey":       "dt_1D",
          "durationLogKey":              "duration_1D",
          "timeStepOutputInterval":      1e4,
          "dirichletBoundaryConditions": bc,
          "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
          "inputMeshIsGlobal":           True,
          "checkForNanInf":              False,
          "solverName":                  "implicitSolver",
          "nAdditionalFieldVariables":   0,
          "additionalSlotNames":         [],
          
          "FiniteElementMethod" : {
            "inputMeshIsGlobal":        True,
            "meshName":                 "MeshFiber"+str(i),
            "prefactor":                Conductivity/(Am*Cm),
            "solverName":               "implicitSolver",
            "slotName":                 "vm",
          },
          "OutputWriter" : [
            {"format": "Paraview",   "outputInterval": output_interval_1D, "filename": variables.outfile_1D.format(i=i), "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering":"incremental"},
            {"format": "PythonFile", "outputInterval": output_interval_1D, "filename": variables.outfile_1D.format(i=i), "binary": True, "onlyNodalValues":False, "fileNumbering":"incremental"}, # also derivatives for hermite
            # {"format": "MegaMol",    "outputInterval": int(1./dt_1D*megamol_output_timestep), "filename": "out/fibers", "timeStepCloseInterval": 7000},
            # {"format": "Paraview",   "outputInterval": output_interval_1D, "filename": variables.outfile_1D.format(i=i)+"_txt", "binary": False, "fixedFormat": False},
            # {"format": "ExFile",     "outputInterval": output_interval_1D, "filename": variables.outfile_1D.format(i=i), "sphereSize": "0.02*0.02*0.02"},
            # {"format": "PythonFile", "outputInterval": output_interval_1D, "filename": variables.outfile_1D.format(i=i), "binary":True, "onlyNodalValues":True},
            # {"format": "PythonFile", "outputInterval": output_interval_1D, "filename": variables.outfile_1D.format(i=i), "binary":False, "onlyNodalValues":True},
            # {"format":"PythonCallback", "outputInterval": output_interval_1D, "callback": lambda x: print("writing fiber {ix:5}: 1D for time step {timeStepNo:5} at time {currentTime}".format(ix=i, **x[i]))},
          ] if variables.outfile_1D != '' else []
        },
      },
    }
  }
  return instance_config
    
    
# create fibre meshes
meshes = {}

# load streamlines from file, either from a .bin file or a python pickle file
if ".bin" in variables.fiber_file:
  # data input from bin files that contain fibers
  try:
    fiber_file_handle = open(variables.fiber_file, "rb")
  except:
    print("Error: Could not open fiber file \"{}\"".format(variables.fiber_file))
    quit()

  # parse fibers from a binary fiber file that was created by parallel_fiber_estimation
  # parse file header to extract number of fibers
  bytes_raw = fiber_file_handle.read(32)
  header_str = struct.unpack('32s', bytes_raw)[0]
  header_length_raw = fiber_file_handle.read(4)
  header_length = struct.unpack('i', header_length_raw)[0]

  parameters = []
  for i in range(int(header_length/4.) - 1):
    double_raw = fiber_file_handle.read(4)
    value = struct.unpack('i', double_raw)[0]
    parameters.append(value)
    
  n_fibers_total = parameters[0]
  n_fibers_x = (int)(np.round(np.sqrt(n_fibers_total)))
  n_fibers_y = n_fibers_x
  n_points_initial_whole_fiber = parameters[1]

  # parse whole fiber file
  streamlines = []
  mesh_node_positions = []
  for fiber_no in range(n_fibers_total):
    fiber = []
    for point_no in range(n_points_initial_whole_fiber):
      point = []
      for i in range(3):
        double_raw = fiber_file_handle.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
      fiber.append(point)
      
    # sample fiber in z direction
    new_fiber = []
    for point_no in range(n_points_initial_whole_fiber):
      point = fiber[point_no]
      new_fiber.append(point)
      mesh_node_positions.append(point)
    
    streamlines.append(new_fiber)
        
else:
  # load pickle file that contains streamlines
  with open(variables.fiber_file, "rb") as f:
    streamlines = pickle.load(f)
    
nInstances = len(streamlines)
if rank_no == 0:
  print("nInstances: {}".format(nInstances))

#nInstances = 1
    
for i,streamline in enumerate(streamlines):
  
  center_node = int(len(streamline)/2)
  #streamline = streamline[center_node-2:center_node+2]
  
  # define mesh
  meshes["MeshFiber{}".format(i)] = {
    "nElements": len(streamline)-1,
    "nodePositions": streamline,
    "inputMeshIsGlobal": True,
    "setHermiteDerivatives": False,
    "logKey": "Fiber{}".format(i)
  }
    
# load MU distribution and firing times
fibre_distribution = np.genfromtxt(variables.fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(variables.firing_times_file)

# determine when the fibres will fire, for debugging output
if rank_no == 0:
  print("Debugging output about fibre firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  
  n_firing_times = np.size(firing_times,0)
  for fibre_no_index in range(nInstances):
    first_stimulation = None
    for current_time in np.linspace(0,1./variables.stimulation_frequency*n_firing_times,n_firing_times):
      if fiber_gets_stimulated(fibre_no_index, variables.stimulation_frequency, current_time):
        first_stimulation = current_time
        break
  
    print("   Fibre {} is of MU {} and will be stimulated for the first time at {}".format(fibre_no_index, get_motor_unit_no(fibre_no_index), first_stimulation))

# create megamol config file
config_file_contents = \
"""print('Hi, I am the megamolconfig.lua!')

-- mmSetAppDir("{megamol_home}/bin")
mmSetAppDir(".")

mmSetLogFile("")
mmSetLogLevel(0)
mmSetEchoLevel('*')

mmAddShaderDir("{megamol_home}/share/shaders")
mmAddResourceDir("{megamol_home}/share/resources")

mmPluginLoaderInfo("{megamol_home}/lib", "*.mmplg", "include")

-- mmSetConfigValue("*-window", "w1280h720")
mmSetConfigValue("*-window", "w720h720")
mmSetConfigValue("consolegui", "on")

mmSetConfigValue("LRHostEnable", "true")

return "done with megamolconfig.lua."
-- error("megamolconfig.lua is not happy!")
""".format(megamol_home="/store/software/opendihu/dependencies/megamol/install")

config_filename = "megamol_config.lua"
with open(config_filename, "w") as f:
  f.write(config_file_contents)

config = {
  #"MegaMolArguments": "--configfile {} -p ../input/testspheres.lua ".format(config_filename),  
  #"MegaMolArguments": "--configfile {} -p ../input/adios_sphere.lua ".format(config_filename),  
  "MegaMolArguments": "--configfile {} -p ../../../input/adios_project.lua ".format(config_filename),  
  #"MegaMolArguments": "--configfile {} -p ../input/empty_view.lua ".format(config_filename),  
  "scenarioName": variables.scenario_name,
  "mappingsBetweenMeshesLogFile": "",
  "logFormat": "csv",
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Meshes": meshes,
  "Solvers": {
    "implicitSolver": {
      "maxIterations":      variables.diffusion_solver_maxit,
      "relativeTolerance":  variables.diffusion_solver_reltol,
      "absoluteTolerance":  0,
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFormat":         "default",
      "dumpFilename": "",
    }
  },
  "MultipleInstances": {
    "nInstances": nInstances,
    "instances": [get_instance_config(i) for i in range(nInstances)],
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fibre_all", "binary": True, "fixedFormat": False, "combineFiles": True},
      #{"format": "MegaMol",  "outputInterval": 1, "filename": "out/fibers", "timeStepCloseInterval": 7000}
      #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
      #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
      #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "binary":True, "onlyNodalValues":True},
    ]
  },
}
