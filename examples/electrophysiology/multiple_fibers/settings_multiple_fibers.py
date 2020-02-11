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

end_time = 50.0

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import sys

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
solver_type = "gmres"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 3e-3                      # overall timestep width of splitting
output_timestep = 4e-1             # timestep for output files
megamol_output_timestep = 5e-2    # timestep frequency for megamol

# input files
#cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../input/shorten.cpp"
cellml_file = "../../input/hodgkin_huxley_1952.c"

#fibre_file = "../../input/laplace3d_structured_quadratic"
fibre_file = "../../input/laplace3d_structured_linear"
#fibre_file = "../../input1000/laplace3d_structured_quadratic"

fibre_distribution_file = "../../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../../input/MU_firing_times_real.txt"
firing_times_file = "../../input/MU_firing_times_immediately.txt"

#print("prefactor: ",Conductivity/(Am*Cm))
#print("numpy path: ",np.__path__)

# parse arguments
if len(sys.argv) == 2:
  n_processes_per_fiber = 1
else:
  try:
    n_processes_per_fiber = (int)(sys.argv[0])
  except:
    n_processes_per_fiber = 1

scenario_name = ""
if len(sys.argv) <= 3:
  scenario_name = ""
else:
  scenario_name = sys.argv[1]


rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

if rank_no == 0:
  print("n_processes_per_fiber: {}".format(n_processes_per_fiber))
  print("scenario_name: {}".format(scenario_name))

#print("rank: {}/{}".format(rank_no,n_ranks))

# set values for cellml model
if "shorten" in cellml_file:
  parameters_used_as_intermediate = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 400.
  
elif "hodgkin_huxley" in cellml_file:
  parameters_used_as_intermediate = []
  parameters_used_as_constant = [2]
  parameters_initial_values = [0.0]
  nodal_stimulation_current = 40.

def get_motor_unit_no(fibre_no):
  return int(fibre_distribution[fibre_no % len(fibre_distribution)]-1)

def fibre_gets_stimulated(fibre_no, frequency, current_time):

  # determine motor unit
  mu_no = (int)(get_motor_unit_no(fibre_no)*0.8)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters_null(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  pass
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(fibre_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)
  
  # stimulation value
  if is_fibre_gets_stimulated:
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
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fibre_no, get_motor_unit_no(fibre_no), is_fibre_gets_stimulated))
    
  #wait = input("Press any key to continue...")
    
# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(fibre_no, stimulation_frequency, current_time)
  
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
  if is_fibre_gets_stimulated:
    stimulation_current = 40.
  else:
    stimulation_current = 0.

  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index)

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(fibre_no, stimulation_frequency, current_time)

  if is_fibre_gets_stimulated:  
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

  bc = {0: -75, -1: -75}
  instance_config = {
    "ranks": ranks,
    "StrangSplitting": {
      #"numberTimeSteps": 1,
      "timeStepWidth": dt_3D,  # 1e-1
      "logTimeStepWidthAsKey": "dt_3D",
      "durationLogKey": "duration_total",
      "timeStepOutputInterval" : 1000,
      "endTime": end_time,
      "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
      "connectedSlotsTerm2To1": [0],   # transfer the same back

      "Term1": {      # CellML
        "Heun" : {
          "timeStepWidth": dt_0D,  # 5e-5
          "logTimeStepWidthAsKey": "dt_0D",
          "durationLogKey": "duration_0D",
          "initialValues": [],
          "timeStepOutputInterval": 1e4,
          "inputMeshIsGlobal": True,
          "dirichletBoundaryConditions": {},
          "nAdditionalFieldVariables": 0,
            
          "CellML" : {
            "modelFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            "compilerFlags":                          "-fPIC -O3 -march=native -shared ",
            "optimizationType":                       "vc",     # "vc", "simd", "openmp" type of generated optimizated source file
            "approximateExponentialFunction":         False,    # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
            "maximumNumberOfThreads":                 0,        # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
            #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
            #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
            "setSpecificStatesCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setSpecificStatesCallFrequency":         0,        # set_specific_states should be called variables.stimulation_frequency times per ms, 0 means disabled
            "setSpecificStatesFrequencyJitter":       0,        # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
            "setSpecificStatesRepeatAfterFirstCall":  0,        # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
            "setSpecificStatesCallEnableBegin":       0,        # [ms] first time when to call setSpecificStates
            
            "additionalArgument": i,
            
            "intermediatesForTransfer":  [],            # which intermediate values to use in further computation
            "statesForTransfer":   0,                   # which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
                      
            "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFiber"+str(i),
            "stimulationLogFilename": "out/stimulation.log",
            
          },
        },
      },
      "Term2": {     # Diffusion
        "ImplicitEuler" : {
          "initialValues": [],
          #"numberTimeSteps": 1,
          "timeStepWidth": dt_1D,  # 1e-5
          "logTimeStepWidthAsKey": "dt_1D",
          "durationLogKey": "duration_1D",
          "timeStepOutputInterval": 1e4,
          "dirichletBoundaryConditions": bc,
          "inputMeshIsGlobal": True,
          "solverName": "implicitSolver",
          "nAdditionalFieldVariables":  0,
          "FiniteElementMethod" : {
            "maxIterations": 1e4,
            "relativeTolerance": 1e-10,
            "inputMeshIsGlobal": True,
            "meshName": "MeshFiber"+str(i),
            "prefactor": Conductivity/(Am*Cm),
            "solverName": "implicitSolver",
          },
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fibre_"+str(i), "binary": True, "fixedFormat": False, "combineFiles": True},
            {"format": "MegaMol",  "outputInterval": int(1./dt_1D*megamol_output_timestep), "filename": "out/fibers", "timeStepCloseInterval": 7000, "combineNInstances": 1, "useFrontBackBuffer": False},
            #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
            #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "sphereSize": "0.02*0.02*0.02"},
            #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "binary":True, "onlyNodalValues":True},
            #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":False, "onlyNodalValues":True},
          ]
        },
      },
    }
  }
  return instance_config
    
    
# create fibre meshes
meshes = {}

with open(fibre_file, "rb") as f:
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
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# determine when the fibres will fire, for debugging output
if rank_no == 0:
  print("Debugging output about fibre firing: Taking input from file \"{}\"".format(firing_times_file))
  
  n_firing_times = np.size(firing_times,0)
  for fibre_no_index in range(nInstances):
    first_stimulation = None
    for current_time in np.linspace(0,1./stimulation_frequency*n_firing_times,n_firing_times):
      if fibre_gets_stimulated(fibre_no_index, stimulation_frequency, current_time):
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
  "MegaMolArguments": "--configfile {} -p ../../input/adios_project.lua ".format(config_filename),  
  #"MegaMolArguments": "--configfile {} -p ../input/empty_view.lua ".format(config_filename),  
  "scenarioName": scenario_name,
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Meshes": meshes,
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": solver_type,
      "preconditionerType": "none",
      "dumpFormat": "default",
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
