# Multicompartment 3D, biceps
#

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

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_0D = 3e-3                      # timestep width of ODEs
output_timestep = 1e0             # timestep for output files
end_time = 1.0                   # end simulation time
#end_time = dt_0D

# input files
#mesh_file = "../input/mesh_tiny"
mesh_file = "../input/mesh_small"
#mesh_file = "../input/mesh_normal"
cellml_file = "../input/hodgkin_huxley_1952.c"
fibre_distribution_file = "../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../input/MU_firing_times_real.txt"
firing_times_file = "../input/MU_firing_times_immediately.txt"

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# load mesh
with open(mesh_file, "rb") as f:
  mesh_data = pickle.load(f, encoding='latin1')
#
#  "node_positions": node_positions, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
#

# load MU distribution and firing times
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# cellml settings
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

def compartmentGetsStimulated(fibre_no, frequency, current_time):
  # determine motor unit
  mu_no = (int)(getMotorUnitNo(fibre_no)*0.8)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, compartment_no):
  
  # determine if fibre gets stimulated at the current time
  compartment_gets_stimulated = compartmentGetsStimulated(compartment_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  
  n_linear_elements_per_coordinate_direction = mesh_data["n_linear_elements_per_coordinate_direction"]
  n_nodes_x = n_linear_elements_per_coordinate_direction[0]+1
  n_nodes_y = n_linear_elements_per_coordinate_direction[1]+1
  n_nodes_z = n_linear_elements_per_coordinate_direction[2]+1
  z_index_center = (int)(n_nodes_z/2)
  
  #nodes_to_stimulate_global = [k*n_nodes_y*n_nodes_x + j*n_nodes_y + i for i in range(n_nodes_x) for j in range(n_nodes_y) for k in [z_index_center-1, z_index_center, z_index_center+1]]

  # stimulation value
  if compartment_gets_stimulated:
    stimulation_current = 400.
  else:
    stimulation_current = 0.
  
  for dof_no_local,dof_no_global in enumerate(dof_nos_global):
    k = (int)(dof_no_global / (n_nodes_x*n_nodes_y))
    j = (int)((dof_no_global % (n_nodes_x*n_nodes_y)) / n_nodes_x)
    i = (int)(dof_no_global % n_nodes_x)
  
    if z_index_center-1 <= k <= z_index_center+1:
      parameters[dof_no_local] = stimulation_current
  
      #print("       {}: set stimulation for local dof {}".format(rank_no, dof_no_local))
  
  #print("       {}: setParameters at timestep {}, t={}, n_nodes_global={}, range: [{},{}], fibre no {}, MU {}, stimulated: {}".\
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fibre_no, getMotorUnitNo(fibre_no), compartment_gets_stimulated))
    
  #wait = input("Press any key to continue...")

# boundary conditions
potential_flow_bc = {}
for bottom_node_index in mesh_data["bottom_nodes"]:
  potential_flow_bc[bottom_node_index] = 0.0
for top_node_index in mesh_data["top_nodes"]:
  potential_flow_bc[top_node_index] = 1.0
  
nCompartments = 1
  
config = {
  "Meshes": {
    "mesh": {
      "nElements": mesh_data["n_linear_elements_per_coordinate_direction"],
      "nodePositions": mesh_data["node_positions"],
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    }
  },
  "Solvers": {
    "potentialFlowSolver": {
      "relativeTolerance": 1e-15,
      "maxIterations": 10000,
      "solverType": "gmres",
      "preconditionerType": "none"
    },
    "activationSolver": {
      "relativeTolerance": 1e-15,
      "maxIterations": 10000,
      "solverType": "gmres",
      "preconditionerType": "none"
    }
  },
  "MultidomainSolver" : {
    "nCompartments": nCompartments,
    "am": Am,
    "cm": Cm,
    "timeStepWidth": dt_0D,
    "endTime": end_time,
    "timeStepOutputInterval": 1,
    "solverName": "activationSolver",
#    "compartmentRelativeFactors": [
#      [...],
#      [...],
#    ],
    "PotentialFlow": {
      "FiniteElementMethod" : {  
        "meshName": "mesh",
        "solverName": "potentialFlowSolver",
        "prefactor": 1.0,
        "dirichletBoundaryConditions": potential_flow_bc,
        "inputMeshIsGlobal": True,
      },
    },
    "Activation": {
      "FiniteElementMethod" : {  
        "meshName": "mesh",
        "solverName": "activationSolver",
        "prefactor": 1.0,
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "diffusionTensor": [                 # fiber direction is (1,0,0)
          1, 0, 0,
          0, 1, 0,
          0, 0, 1
        ],
        "extracellularDiffusionTensor": [
          2, 0, 0,
          0, 1, 0,
          0, 0, 1
        ]
      },
    },
    "CellMLAdapters": [
      {
        "CellML" : {
          "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
          #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
          #"libraryFilename": "cellml_simd_lib.so",   # compiled library
          "useGivenLibrary": False,
          #"statesInitialValues": [],
          "setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
          "setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setParametersFunctionAdditionalParameter": compartment_no,
          
          "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
          "parametersUsedAsIntermediate": parametersUsedAsIntermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersUsedAsConstant": parametersUsedAsConstant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersInitialValues": parametersInitialValues,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "meshName": "mesh",
          "prefactor": 1.0,
        }
      } for compartment_no in range(nCompartments)
    ],
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/output", "binary": True, "fixedFormat": False, "combineFiles": False},
      #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
      #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
    ]
  },
}
