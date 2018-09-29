# Cuboid
#
# arguments: <n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name>

end_time = 10.0

import numpy as np
import pickle
import sys

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
  
cellml_file = "../input/shorten_ocallaghan_davidson_soboleva_2007.c"
cellml_file = "../input/shorten_opencmiss.cpp"
cellml_file = "../input/hodgkin_huxley_1952.c"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 3e-3                      # overall timestep width of splitting
output_timestep = 1e0             # timestep for output files
n_nodes_per_fiber = 1000             # number of nodes per fiber
n_fibers = 10

#print("prefactor: ",Conductivity/(Am*Cm))
#print("numpy path: ",np.__path__)

# parse arguments
n_processes_per_fiber = (int)(sys.argv[0])
n_fibers = (int)(sys.argv[1])
n_nodes_per_fiber = (int)(sys.argv[2])
scenario_name = sys.argv[3]

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

if rank_no == 0:
  print("scenario_name: {}".format(scenario_name))
  print("n_processes_per_fiber: {}, n_fibers: {}, n_nodes_per_fiber: {}".format(n_processes_per_fiber, n_fibers, n_nodes_per_fiber))

#print("rank: {}/{}".format(rank_no,n_ranks))

if "shorten" in cellml_file:
  parametersUsedAsIntermediate = [32]
  parametersUsedAsConstant = [65]
  parametersInitialValues = [0.0, 1.0]
  
elif "hodgkin_huxley" in cellml_file:
  parametersUsedAsIntermediate = []
  parametersUsedAsConstant = [2]
  parametersInitialValues = [0.0]
  
def fibreGetsStimulated(fibre_no, frequency, current_time):

  # determine if fibre fires now
  index = int(current_time * frequency)
  if index % 10 == 0:
    return True
  else:
    return False
  
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
  if n_processes_per_fiber > 0:
    ranks = []
    for j in range(n_processes_per_fiber):
      ranks.append(n_processes_per_fiber*i + j)
  else:
    ranks = [int(i/-n_processes_per_fiber)]

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
      "outputData1": False,
      "outputData2": True,

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
            "meshName": "MeshFiber",
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
            "meshName": "MeshFiber",
            "prefactor": Conductivity/(Am*Cm),
          },
          "OutputWriter" : [
            #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i), "binary": True, "fixedFormat": False},
            #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
            #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
            #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
          ]
        },
      },
    }
  }
  return instance_config
    
config = {
  "scenarioName": scenario_name,
  "Meshes": {
    "MeshFiber": {
      "nElements": n_nodes_per_fiber-1,
      "nodePositions": [[x,0,0] for x in np.linspace(0,(n_nodes_per_fiber-1)/100.,n_nodes_per_fiber)],
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False,
      "logKey": "1D"
    }
  },
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
    }
  },
  "MultipleInstances": {
    "nInstances": n_fibers,
    "instances": [get_instance_config(i) for i in range(n_fibers)],
  }
}
