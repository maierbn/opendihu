# Multicompartment 3D, biceps
#

import numpy as np
import scipy.stats
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
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 3e-3                      # overall timestep width of splitting
output_timestep = 1e-2             # timestep for output files
end_time = 500.0                   # end simulation time
#end_time = dt_0D

Am = 0.2   # mesh_small

# input files
#mesh_file = "../input/scaled_mesh_tiny"
mesh_file = "../input/scaled_mesh_small"
#mesh_file = "../input/scaled_mesh_normal"
#mesh_file = "../input/scaled_mesh_big"
fiber_file = "../input/laplace3d_structured_linear"
cellml_file = "../input/hodgkin_huxley_1952.c"
fibre_distribution_file = "../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../input/MU_firing_times_real.txt"
firing_times_file = "../input/MU_firing_times_immediately.txt"

# motor unit parameters
motor_units = [
  {"fiber_no": 10, "standard_deviation": 2.0, "maximum": 0.5},
  {"fiber_no": 30, "standard_deviation": 2.0, "maximum": 0.4},
  {"fiber_no": 50, "standard_deviation": 3.0, "maximum": 0.6},
]


if "hodgkin" in cellml_file:
  Cm = 1.0

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

# load fibers
with open(fiber_file, "rb") as f:
  fiber_data = pickle.load(f, encoding='latin1')
# list of fibers, fiber = list of points, point = list with 3 coordinate entries

min_x = min([x for [x,y,z] in mesh_data["node_positions"]])
max_x = max([x for [x,y,z] in mesh_data["node_positions"]])
min_y = min([y for [x,y,z] in mesh_data["node_positions"]])
max_y = max([y for [x,y,z] in mesh_data["node_positions"]])
min_z = min([z for [x,y,z] in mesh_data["node_positions"]])
max_z = max([z for [x,y,z] in mesh_data["node_positions"]])

if rank_no == 0:
  print("mesh bounding box x: [{},{}], y: [{},{}], z:[{},{}]".format(min_x, max_x, min_y, max_y, min_z, max_z))

for fiber_no in [10, 30, 50]:
  data = fiber_data[fiber_no]
  min_x = min([x for [x,y,z] in data])
  max_x = max([x for [x,y,z] in data])
  min_y = min([y for [x,y,z] in data])
  max_y = max([y for [x,y,z] in data])
  min_z = min([z for [x,y,z] in data])
  max_z = max([z for [x,y,z] in data])

  if rank_no == 0:
    print("fiber {} bounding box x: [{},{}], y: [{},{}], z:[{},{}]".format(fiber_no, min_x, max_x, min_y, max_y, min_z, max_z))

n_compartments = len(motor_units)

# create relative factors for compartments

if rank_no == 0:
  print("determine relative factors for {} motor units:\n{}".format(n_compartments, motor_units))

# create data structure with 0
relative_factors = np.zeros((n_compartments, len(mesh_data["node_positions"])))   # each row is one compartment

# loop over nodes of mesh
for node_no,node_position in enumerate(mesh_data["node_positions"]):
  node_position = np.array(node_position)
  
  # loop over motor units
  for motor_unit_no,motor_unit in enumerate(motor_units):
    
    # find point on fiber that is closest to current node
    fiber_no = motor_unit["fiber_no"]
    if fiber_no >= len(fiber_data):
      print("Error with motor unit {}, only {} fibers available".format(motor_unit, len(fiber_datar)))
    else:
      max_distance = None
      for fiber_point in fiber_data[fiber_no]:
        d = np.array(fiber_point) - node_position
        distance = np.inner(d,d)
        if max_distance is None or distance < max_distance:
          max_distance = distance
          #print("node_position {}, fiber_point {}, d={}, |d|={}".format(node_position, fiber_point, d, np.sqrt(distance)))
      
      distance = np.sqrt(max_distance)
      
      
      gaussian = scipy.stats.norm(loc = 0., scale = motor_unit["standard_deviation"])
      value = gaussian.pdf(distance)*motor_unit["maximum"]
      relative_factors[motor_unit_no][node_no] += value
      #print("motor unit {}, fiber {}, distance {}, value {}".format(motor_unit_no, fiber_no, distance, value))

if rank_no == 0:
  for i,factors_list in enumerate(relative_factors.tolist()):
    print("MU {}, maximum fr: {}".format(i,max(factors_list)))

# load MU distribution and firing times
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# cellml settings
if "shorten" in cellml_file:
  parameters_used_as_intermediate = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  
elif "hodgkin_huxley" in cellml_file:
  parameters_used_as_intermediate = []
  parameters_used_as_constant = [2]
  parameters_initial_values = [0.0]
  
def get_motor_unit_no(fibre_no):
  return int(fibre_distribution[fibre_no % len(fibre_distribution)]-1)

def compartment_gets_stimulated(compartment_no, current_time):
  # determine motor unit
  mu_no = (int)(get_motor_unit_no(compartment_no)*0.8)
  
  # determine if MU fires now
  index = int(current_time * stimulation_frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, compartment_no):
  
  # determine if fibre gets stimulated at the current time
  is_compartment_gets_stimulated = compartment_gets_stimulated(compartment_no, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  
  n_linear_elements_per_coordinate_direction = mesh_data["n_linear_elements_per_coordinate_direction"]
  n_nodes_x = n_linear_elements_per_coordinate_direction[0]+1
  n_nodes_y = n_linear_elements_per_coordinate_direction[1]+1
  n_nodes_z = n_linear_elements_per_coordinate_direction[2]+1
  z_index_center = (int)(n_nodes_z/2)
  y_index_center = (int)(n_nodes_y/2)
  x_index_center = (int)(n_nodes_x/2)
  
  #nodes_to_stimulate_global = [k*n_nodes_y*n_nodes_x + j*n_nodes_y + i for i in range(n_nodes_x) for j in range(n_nodes_y) for k in [z_index_center-1, z_index_center, z_index_center+1]]

  # stimulation value
  if is_compartment_gets_stimulated:
    stimulation_current = 400.
  else:
    stimulation_current = 0.
  
  for dof_no_local,dof_no_global in enumerate(dof_nos_global):
    k = (int)(dof_no_global / (n_nodes_x*n_nodes_y))
    j = (int)((dof_no_global % (n_nodes_x*n_nodes_y)) / n_nodes_x)
    i = (int)(dof_no_global % n_nodes_x)
  
    if z_index_center-1 <= k <= z_index_center+1:
      if y_index_center-1 <= j <= y_index_center+1:
        if x_index_center-1 <= i <= x_index_center+1:
      
          parameters[dof_no_local] = stimulation_current
  
      #print("       {}: set stimulation for local dof {}".format(rank_no, dof_no_local))
  
  #print("       {}: setParameters at timestep {}, t={}, n_nodes_global={}, range: [{},{}], fibre no {}, MU {}, stimulated: {}".\
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fibre_no, get_motor_unit_no(fibre_no), compartment_gets_stimulated))
    
  #wait = input("Press any key to continue...")

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, compartment_no):
  
  # determine if fibre gets stimulated at the current time
  is_compartment_gets_stimulated = compartment_gets_stimulated(compartment_no, current_time)

  if is_compartment_gets_stimulated:  
    
    n_linear_elements_per_coordinate_direction = mesh_data["n_linear_elements_per_coordinate_direction"]
    n_nodes_x = n_linear_elements_per_coordinate_direction[0]+1
    n_nodes_y = n_linear_elements_per_coordinate_direction[1]+1
    n_nodes_z = n_linear_elements_per_coordinate_direction[2]+1
    z_index_center = (int)(n_nodes_z/2)
    y_index_center = (int)(n_nodes_y/2)
    x_index_center = (int)(n_nodes_x/2)
    
    for k in range(n_nodes_z):
      if z_index_center-1 <= k <= z_index_center+1:
        for j in range(n_nodes_y):
          if y_index_center-1 <= j <= y_index_center+1:
            for i in range(n_nodes_x):
              if x_index_center-1 <= i <= x_index_center+1:
                key = ((i,j,k),0,0)        # key: ((x,y,z),nodal_dof_index,state_no)
                states[key] = 400.0
                #print("set states at ({},{},{}) to 400".format(i,j,k))

    #print("states: {}".format(states))
    #print("n_nodes: ({},{},{})".format(n_nodes_x, n_nodes_y, n_nodes_z))
    #print("n_nodes_global: {}, time_step_no: {}, current_time: {}, compartment_no: {}".format(n_nodes_global, time_step_no, current_time, compartment_no))
    #wait = input("Press any key to continue...")
    
# boundary conditions
potential_flow_bc = {}
for bottom_node_index in mesh_data["bottom_nodes"]:
  potential_flow_bc[bottom_node_index] = 0.0
for top_node_index in mesh_data["top_nodes"]:
  potential_flow_bc[top_node_index] = 1.0
  
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
      "relativeTolerance": 1e-10,
      "maxIterations": 10000,
      "solverType": "gmres",
      "preconditionerType": "none"
    },
    "activationSolver": {
      "relativeTolerance": 1e-5,
      "maxIterations": 10000,
      "solverType": "gmres",
      "preconditionerType": "none"
    }
  },
  "StrangSplitting": {
    "timeStepWidth": dt_3D,  # 1e-1
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval" : 100,
    "endTime": end_time,
    "Term1": {      # CellML
      "MultipleInstances": {
        "nInstances": n_compartments,  
        "instances": [        # settings for each motor unit, `i` is the index of the motor unit
        {
          "ranks": list(range(n_ranks)),
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
              "compilerFlags": "-fPIC -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared ",
              "useGivenLibrary": False,
              #"statesInitialValues": [],
              "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              "setSpecificStatesCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
              "additionalArgument": compartment_no,
              #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
              #"setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
              #"setParametersFunctionAdditionalParameter": compartment_no,
              
              "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
              "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
              "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
              "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
              "meshName": "mesh",
              "prefactor": 1.0,
            }
          }
        } for compartment_no in range(n_compartments)]
      },
    },
    "Term2": {     # Diffusion
      "MultidomainSolver" : {
        "nCompartments": n_compartments,
        "am": Am,
        "cm": Cm,
        "timeStepWidth": dt_0D,
        "endTime": end_time,
        "timeStepOutputInterval": 50,
        "solverName": "activationSolver",
        "inputIsGlobal": True,
        "compartmentRelativeFactors": relative_factors.tolist(),
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
        
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": (int)(1./dt_1D*output_timestep), "filename": "out/output", "binary": True, "fixedFormat": False, "combineFiles": False},
          #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
          #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
        ]
      }
    }
  }
}
