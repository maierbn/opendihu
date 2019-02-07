# 2 fibers, biceps, example for load_balancing
#

end_time = 150.0     # end time for the simulation

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
solver_type = "gmres"   # solver for the linear system

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 0.1                         # overall timestep width of splitting
output_timestep = 1e0            # timestep for output files

# input files
fibre_file = "../input/laplace3d_structured_linear"
fibre_distribution_file = "../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../input/MU_firing_times_load_balancing.txt"
cellml_file = "../input/hodgkin_huxley_1952.c"

# get own rank no and number of ranks
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

def getMotorUnitNo(fibre_no):
  """
  get the no. of the motor unit which fiber fibre_no is part of
  """
  return (int)(0.8*int(fibre_distribution[fibre_no % len(fibre_distribution)]-1))

def fibreGetsStimulated(fibre_no, frequency, current_time):
  """
  determine if fiber fibre_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  mu_no = getMotorUnitNo(fibre_no)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  """
  This function sets the stimulation current when the fiber will fire.
  :param n_nodes_global: number of global nodes of the fiber
  :param time_step_no:  the index of the time step in the simulation, integer value
  :param current_time:  the current simulation time as double
  :param parameters:    a list of parameter values for the nodes, changes made to the list will be used by the simulation
  :param dof_nos_global: for every local degree of freedom no. the global dof no., i.e. this is a local to global mapping
  :param fibre_no:      the no. of the fiber (here 0 or 1). It is the value that was set by "setParametersFunctionAdditionalParameter"
  """ 
  
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
    
# load fiber meshes, called streamlines
with open(fibre_file, "rb") as f:
  streamlines = pickle.load(f)

# assign fiber meshes to names streamline0, streamline1
if len(streamlines) < 2:
  print("Error: input file {} only contains {} fibers, 2 needed".format(fibre_file, len(streamlines)))
streamline0 = streamlines[0]
streamline1 = streamlines[1]
    
# load MU distribution and firing times
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# determine when the fibres will fire, this is only for debugging output and is not used by the config
if rank_no == 0:
  print("Debugging output about fibre firing: Taking input from file \"{}\"".format(firing_times_file))
  
  n_firing_times = np.size(firing_times,0)
  for fibre_no_index in range(2):
    first_stimulation = None
    for current_time in np.linspace(0,1./stimulation_frequency*n_firing_times,n_firing_times):
      if fibreGetsStimulated(fibre_no_index, stimulation_frequency, current_time):
        first_stimulation = current_time
        break
  
    mu_no = getMotorUnitNo(fibre_no_index)
    print("   Fibre {} is of MU {} and will be stimulated for the first time at {}".format(fibre_no_index, mu_no, first_stimulation))

# configure on which ranks fibers 0 and 1 will run
ranks = [
  [0,1],       # rank nos that will compute fiber 0
  [2,3]        # rank nos that will compute fiber 1
]

config = {
  "scenarioName": "load_balancing",
  # the 1D meshes for each fiber
  "Meshes": {
    "MeshFibre0":
    {
      "nElements": len(streamline0)-1,
      "nodePositions": streamline0,
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    },
    "MeshFibre1":
    {
      "nElements": len(streamline1)-1,
      "nodePositions": streamline1,
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    }
  },
  # the linear solver for the diffusion problem
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": solver_type,
      "preconditionerType": "none",
    }
  },
  # control class that configures multiple instances of the fiber model
  "MultipleInstances": {
    "nInstances": 2,      # number of fibers
    "instances": [        # settings for each fiber, `i` is the index of the fiber (0 or 1)
    {
      "ranks": ranks[i],
      
      # config for strang splitting
      "StrangSplitting": {
        "timeStepWidth": dt_3D,  # 1e-1
        "logTimeStepWidthAsKey": "dt_3D",
        "durationLogKey": "duration_total",
        "timeStepOutputInterval" : 1,
        "endTime": end_time,

        "Term1": {      # CellML
          "LoadBalancing": {
            "Heun" : {
              "timeStepWidth": dt_0D,  # 5e-5
              "maxTimeStepWidth": 0.1,
              "tolerance": 0.001,
              "logTimeStepWidthAsKey": "dt_0D",
              "durationLogKey": "duration_0D",
              "initialValues": [],
              "timeStepOutputInterval": 1,
              "inputMeshIsGlobal": True,
              "dirichletBoundaryConditions": {},
                
              "CellML" : {
                "sourceFilename": cellml_file,                     # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
                "useGivenLibrary": False,
                "setParametersFunction": set_parameters,           # callback function that sets parameters like stimulation current
                "setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # interval in which the set_parameters function will be called
                "setParametersFunctionAdditionalParameter": i,
                
                "outputStateIndex": 0,                             # state 0 = Vm, rate 28 = gamma
                "parametersUsedAsIntermediate": [],                # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                "parametersUsedAsConstant": [2],                   # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                "parametersInitialValues": [0.0],                  # initial values for the parameters
                "meshName": "MeshFibre"+str(i),                    # name of the fiber mesh, i.e. either MeshFibre0 or MeshFibre1
                "prefactor": 1.0,
              },
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
            "dirichletBoundaryConditions": {0: -75, -1: -75},      # set first and last value of fiber to -75
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
    for i in range(2)
    ]
  }
}
