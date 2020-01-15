# Multicompartment 3D, biceps
#

import numpy as np
import scipy.stats
import pickle
import sys,os 

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
output_timestep = 1e-1             # timestep for output files
end_time = 500.0                   # end simulation time
#end_time = dt_0D

Am = 0.2   # mesh_small
Am = 0.1

# input files
#mesh_file = "../input/scaled_mesh_tiny"
#mesh_file = "../input/scaled_mesh_small"
mesh_file = "../input/scaled_mesh_normal"
#mesh_file = "../input/scaled_mesh_big"
fiber_file = "../input/laplace3d_structured_linear"

#fiber_file = "../../input/7x7fibers.bin"

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

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])


if ".bin" in fiber_file:
  # data input from bin files that contain fibers

  try:
    fiber_file_handle = open(fiber_file, "rb")
  except:
    print("Error: Could not open fiber file \"{}\"".format(fiber_file))
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
  n_points_whole_fiber = parameters[1]

  n_linear_elements_per_coordinate_direction = [n_fibers_x-1, n_fibers_y-1, n_points_whole_fiber-1]

  if rank_no == 0:
    print("n fibers:              {} ({} x {})".format(n_fibers_total, n_fibers_x, n_fibers_y))
    print("n points per fiber:    {}".format(n_points_whole_fiber))
    
  # parse whole fiber file
  fiber_data = []
  mesh_node_positions = []
  for fiber_no in range(n_fibers_total):
    fiber = []
    for point_no in range(n_points_whole_fiber):
      point = []
      for i in range(3):
        double_raw = fiber_file_handle.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
      fiber.append(point)
      mesh_node_positions.append(point)
    fiber_data.append(fiber)
            
  bottom_node_indices = list(range(n_fibers_x*n_fibers_y))
  n_points = n_fibers_x*n_fibers_y*n_points_whole_fiber
  top_node_indices = list(range(n_points-n_fibers_x*n_fibers_y,n_points))

else:
  # data input from generating 3D meshes without fiber tracing  
  # load fibers
  with open(fiber_file, "rb") as f:
    fiber_data = pickle.load(f, encoding='latin1')
  # list of fibers, fiber = list of points, point = list with 3 coordinate entries

  # load mesh
  with open(mesh_file, "rb") as f:
    mesh_data = pickle.load(f, encoding='latin1')

  n_linear_elements_per_coordinate_direction = mesh_data["n_linear_elements_per_coordinate_direction"]
  mesh_node_positions = mesh_data["node_positions"]

  bottom_node_indices = mesh_data["bottom_nodes"]
  top_node_indices = mesh_data["top_nodes"]

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

  # output bounding box for debugging
  if rank_no == 0:
    min_x = min([x for [x,y,z] in mesh_data["node_positions"]])
    max_x = max([x for [x,y,z] in mesh_data["node_positions"]])
    min_y = min([y for [x,y,z] in mesh_data["node_positions"]])
    max_y = max([y for [x,y,z] in mesh_data["node_positions"]])
    min_z = min([z for [x,y,z] in mesh_data["node_positions"]])
    max_z = max([z for [x,y,z] in mesh_data["node_positions"]])

    print("mesh bounding box x: [{},{}], y: [{},{}], z:[{},{}]".format(min_x, max_x, min_y, max_y, min_z, max_z))

    for fiber_no in [10, 30, 50]:
      data = fiber_data[fiber_no]
      min_x = min([x for [x,y,z] in data])
      max_x = max([x for [x,y,z] in data])
      min_y = min([y for [x,y,z] in data])
      max_y = max([y for [x,y,z] in data])
      min_z = min([z for [x,y,z] in data])
      max_z = max([z for [x,y,z] in data])

      print("fiber {} bounding box x: [{},{}], y: [{},{}], z:[{},{}]".format(fiber_no, min_x, max_x, min_y, max_y, min_z, max_z))

# determine relative factor fields fr(x) for compartments
relative_factors_file = "{}.compartment_relative_factors".format(mesh_file)
if os.path.exists(relative_factors_file):
  with open(relative_factors_file, "rb") as f:
    if rank_no == 0:
      print("load relative factors from file \"{}\"".format(relative_factors_file))
    relative_factors = pickle.load(f, encoding='latin1')

else:
  sys.path.append(os.path.abspath(".."))
  import initialize_compartment_relative_factors
  relative_factors = initialize_compartment_relative_factors.compute_compartment_relative_factors(mesh_node_positions, fiber_data, motor_units)
  if rank_no == 0:
    print("save relative factors to file \"{}\"".format(relative_factors_file))
    with open(relative_factors_file, "wb") as f:
      pickle.dump(relative_factors, f)
  
if "hodgkin" in cellml_file:
  Cm = 1.0

n_compartments = len(motor_units)

# create relative factors for compartments

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
for bottom_node_index in bottom_node_indices:
  potential_flow_bc[bottom_node_index] = 0.0
for top_node_index in top_node_indices:
  potential_flow_bc[top_node_index] = 1.0
  
multidomain_solver = {
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
      "neumannBoundaryConditions": [],
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
      "neumannBoundaryConditions": [],
      "diffusionTensor": [      # sigma_i           # fiber direction is (1,0,0)
        8.93, 0, 0,
        0, 0.893, 0,
        0, 0, 0.893
      ], 
      "extracellularDiffusionTensor": [      # sigma_e
        6.7, 0, 0,
        0, 6.7, 0,
        0, 0, 6.7,
      ],
    },
  },
  
  "OutputWriter" : [
    {"format": "Paraview", "outputInterval": (int)(1./dt_1D*output_timestep), "filename": "out/output", "binary": True, "fixedFormat": False, "combineFiles": False},
    #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
    #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
  ]
}
  
config = {
  "Meshes": {
    "mesh": {
      "nElements": n_linear_elements_per_coordinate_direction,
      "nodePositions": mesh_node_positions,
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    }
  },
  "Solvers": {
    "potentialFlowSolver": {
      "relativeTolerance": 1e-10,
      "maxIterations": 1e4,
      "solverType": "gmres",
      "preconditionerType": "none",
      "dumpFormat": "default",
      "dumpFilename": "",
    },
    "activationSolver": {
      "relativeTolerance": 1e-5,
      "maxIterations": 1e5,
      "solverType": "gmres",
      "preconditionerType": "none",
      "dumpFormat": "default",
      "dumpFilename": "",
    }
  },
  "StrangSplitting": {
    "timeStepWidth": dt_3D,  # 1e-1
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval" : 10,
    "endTime": end_time,
    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
    "connectedSlotsTerm2To1": [0],   # transfer the same back

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
            "nAdditionalFieldVariables": 0,
                
            "CellML" : {
              "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
              #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
              #"libraryFilename": "cellml_simd_lib.so",   # compiled library
              "compilerFlags": "-fPIC -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared ",
              "useGivenLibrary": False,
              #"statesInitialValues": [],
              "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              "setSpecificStatesCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
              "setSpecificStatesCallFrequency":     0,          # set_specific_states should be called variables.stimulation_frequency times per ms
              "setSpecificStatesFrequencyJitter":   0,          # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
              "setSpecificStatesRepeatAfterFirstCall":  0.01,   # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
              "setSpecificStatesCallEnableBegin":   0,          # [ms] first time when to call setSpecificStates
              "additionalArgument": compartment_no,
              #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
              #"setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
              #"setParametersFunctionAdditionalParameter": compartment_no,
              "intermediatesForTransfer":  [],            # which intermediate values to use in further computation
              "statesForTransfer": 0,                     # which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
                     
              "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
              "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
              "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
              "meshName": "mesh",
              "stimulationLogFilename": "out/stimulation.log",
            }
          }
        } for compartment_no in range(n_compartments)]
      },
    },
    "Term2": {     # Diffusion
      "MultidomainSolver" : multidomain_solver,
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": (int)(1./dt_1D*output_timestep), "filename": "out/surface", "binary": True, "fixedFormat": False, "combineFiles": True},
        ],
        "face": "1-",
        "MultidomainSolver" : multidomain_solver,
      }
    }
  }
}
