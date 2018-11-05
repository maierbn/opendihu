# multiple fibers, biceps
#

end_time = 10.0

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
  
cellml_file = "../input/hodgkin_huxley_1952.c"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 3e-3                      # overall timestep width of splitting
output_timestep = 1e0             # timestep for output files

# input files
#fibre_file = "../input/laplace3d_structured_quadratic"
fibre_file = "../input/laplace3d_structured_linear"
#fibre_file = "../input1000/laplace3d_structured_quadratic"
fibre_distribution_file = "../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../input/MU_firing_times_real.txt"
firing_times_file = "../input/MU_firing_times_immediately.txt"

#print("prefactor: ",Conductivity/(Am*Cm))
#print("numpy path: ",np.__path__)

# parse arguments
if len(sys.argv) == 2:
  n_processes_per_fiber = 1
else:
  n_processes_per_fiber = (int)(sys.argv[0])

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

if rank_no == 0:
  print("n_processes_per_fiber: {}".format(n_processes_per_fiber))

#print("rank: {}/{}".format(rank_no,n_ranks))

if "shorten" in cellml_file:
  parametersUsedAsIntermediate = [32]
  parametersUsedAsConstant = [65]
  parametersInitialValues = [0.0, 1.0]
  
elif "hodgkin_huxley" in cellml_file:
  parametersUsedAsIntermediate = []
  parametersUsedAsConstant = [2]
  parametersInitialValues = [0.0]
  

def getMotorUnitNo(fibre_no):
  return int(fibre_distribution[fibre_no % len(fibre_distribution)]-1)

def fibreGetsStimulated(fibre_no, frequency, current_time):

  # determine motor unit
  mu_no = (int)(getMotorUnitNo(fibre_no)*0.8)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters_null(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  pass
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  fibre_gets_stimulated = fibreGetsStimulated(fibre_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)
  
  # stimulation value
  if fibre_gets_stimulated:
    stimulation_current = 400.
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
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fibre_no, getMotorUnitNo(fibre_no), fibre_gets_stimulated))
    
  #wait = input("Press any key to continue...")
    
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

      "Term1": {      # CellML
        "Heun" : {
          "timeStepWidth": dt_0D,  # 5e-5
          "logTimeStepWidthAsKey": "dt_0D",
          "durationLogKey": "duration_0D",
          "initialValues": [],
          "timeStepOutputInterval": 1e4,
          "inputMeshIsGlobal": True,
          "dirichletBoundaryConditions": {},
            
          "CellML" : {
            "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
            #"libraryFilename": "cellml_simd_lib.so",   # compiled library
            "useGivenLibrary": False,
            #"statesInitialValues": [],
            "setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
            "setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setParametersFunctionAdditionalParameter": i,
            
            "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
            "parametersUsedAsIntermediate": parametersUsedAsIntermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parametersUsedAsConstant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parametersInitialValues,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFibre"+str(i),
            "prefactor": 1.0,
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
          "FiniteElementMethod" : {
            "maxIterations": 1e4,
            "relativeTolerance": 1e-10,
            "inputMeshIsGlobal": True,
            "meshName": "MeshFibre"+str(i),
            "prefactor": Conductivity/(Am*Cm),
            "solverName": "implicitSolver",
          },
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fibre_"+str(i), "binary": True, "fixedFormat": False, "combineFiles": False},
            #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
            #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
            {"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "binary":True, "onlyNodalValues":True},
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
  meshes["MeshFibre{}".format(i)] = {
    "nElements": len(streamline)-1,
    "nodePositions": streamline,
    "inputMeshIsGlobal": True,
    "setHermiteDerivatives": False
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
      if fibreGetsStimulated(fibre_no_index, stimulation_frequency, current_time):
        first_stimulation = current_time
        break
  
    print("   Fibre {} is of MU {} and will be stimulated for the first time at {}".format(fibre_no_index, getMotorUnitNo(fibre_no_index), first_stimulation))

config = {
  "scenarioName": "iobuf",
  "Meshes": meshes,
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": solver_type,
      "preconditionerType": "none",
    }
  },
  "MultipleInstances": {
    "nInstances": nInstances,
    "instances": [get_instance_config(i) for i in range(nInstances)],
  }
}
