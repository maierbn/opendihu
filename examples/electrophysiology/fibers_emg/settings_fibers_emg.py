# multiple fibers, biceps
#

end_time = 10.0

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import sys
import struct

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
dt_bidomain = 1e-3                # time step width for bidomain equation
output_timestep = 4e-1             # timestep for output files

# input files
#cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../input/shorten.cpp"
cellml_file = "../../input/hodgkin_huxley_1952.c"

fibre_file = "../../input/3000fibers.bin"
fibre_file = "../../input/500fibers.bin"
#fibre_file = "../../input/49fibers.bin"

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
    
# create fibre meshes
meshes = {}

# parse fibers from a binary fiber file that was created by parallel_fiber_estimation
with open(fibre_file, "rb") as f:

  # parse header
  bytes_raw = f.read(32)
  header_str = struct.unpack('32s', bytes_raw)[0]
  header_length_raw = f.read(4)
  header_length = struct.unpack('i', header_length_raw)[0]
  
  parameters = []
  for i in range(int(header_length/4.) - 1):
    double_raw = f.read(4)
    value = struct.unpack('i', double_raw)[0]
    parameters.append(value)
    
  n_fibers_total = parameters[0]
  n_fibers_x = (int)(np.round(np.sqrt(n_fibers_total)))
  n_fibers_y = n_fibers_x
  n_points_whole_fiber = parameters[1]
  
  # for debugging:
  #n_fibers_total = 4
  #n_fibers_x = 2
  #n_fibers_y = 2
  
  if rank_no == 0:
    print("nFibersTotal:      {} ({} x {})".format(n_fibers_total, n_fibers_x, n_fibers_y))
    print("nPointsWholeFiber: {}".format(n_points_whole_fiber))
  
  fibers = []
  for fiber_no in range(n_fibers_total):
    fiber = []
    for point_no in range(n_points_whole_fiber):
      point = []
      for i in range(3):
        double_raw = f.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
      fiber.append(point)
    fibers.append(fiber)
      
# determine number of fibers and create meshes
if rank_no == 0:
  print("n_fibers_total: {}".format(n_fibers_total))
    
for i,fiber in enumerate(fibers):
  
  center_node = int(len(fiber)/2)
  #fiber = fiber[center_node-2:center_node+2]
  
  # define mesh
  meshes["MeshFiber_{}".format(i)] = {
    "nElements": len(fiber)-1,
    "nodePositions": fiber,
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
  for fibre_no_index in range(n_fibers_total):
    first_stimulation = None
    for current_time in np.linspace(0,1./stimulation_frequency*n_firing_times,n_firing_times):
      if fibre_gets_stimulated(fibre_no_index, stimulation_frequency, current_time):
        first_stimulation = current_time
        break
  
    print("   Fibre {} is of MU {} and will be stimulated for the first time at {}".format(fibre_no_index, get_motor_unit_no(fibre_no_index), first_stimulation))

# compute partitioning
# this has to match the total number of processes
n_subdomains_x = 2   # example values for 4 processes
n_subdomains_y = 1
n_subdomains_z = 2

# here any number is possible
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 3

if rank_no == 0:
  if n_ranks != n_subdomains_x*n_subdomains_y*n_subdomains_z:
    print("Error! Number of ranks {} does not match given partitioning {} x {} x {} ".format(n_ranks, n_subdomains_x, n_subdomains_y, n_subdomains_z))
  
n_subdomains_xy = n_subdomains_x * n_subdomains_y
n_fibers_per_subdomain_x = (int)(np.ceil(n_fibers_x / n_subdomains_x))
n_fibers_per_subdomain_y = (int)(np.ceil(n_fibers_y / n_subdomains_y))
n_points_per_subdomain_z = (int)(np.ceil(n_points_whole_fiber / n_subdomains_z))

if rank_no == 0:
  print("{} ranks, partitioning: x{} x y{} x z{}".format(n_ranks, n_subdomains_x, n_subdomains_y, n_subdomains_z))
  print("{} x {} fibers, per partition: {} x {}, {} points per fiber".format(n_fibers_x, n_fibers_y, n_fibers_per_subdomain_x, n_fibers_per_subdomain_y, n_points_whole_fiber))

# define helper functions for fiber numbering

# number of fibers that are handled inside the subdomain x
def n_fibers_in_subdomain_x(subdomain_coordinate_x):
  return min(n_fibers_per_subdomain_x, n_fibers_x-subdomain_coordinate_x*n_fibers_per_subdomain_x)
  
# number of fibers that are handled inside the subdomain y
def n_fibers_in_subdomain_y(subdomain_coordinate_y):
  return min(n_fibers_per_subdomain_y, n_fibers_y-subdomain_coordinate_y*n_fibers_per_subdomain_y)

# global fiber no, from subdomain coordinate and coordinate inside the subdomain
def fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y):
  return (subdomain_coordinate_y*n_fibers_per_subdomain_y + fiber_in_subdomain_coordinate_y)*n_fibers_x + subdomain_coordinate_x*n_fibers_per_subdomain_x + fiber_in_subdomain_coordinate_x

# number of points that are handled inside the subdomain z
def n_points_in_subdomain_z(subdomain_coordinate_z):
  return min(n_points_per_subdomain_z, n_points_whole_fiber-subdomain_coordinate_z*n_points_per_subdomain_z)
  
def n_sampled_points_in_subdomain_x(subdomain_coordinate_x):
  result = (int)(np.ceil(n_fibers_in_subdomain_x(subdomain_coordinate_x) / sampling_stride_x))
  if subdomain_coordinate_x == n_subdomains_x-1 and (n_fibers_in_subdomain_x(subdomain_coordinate_x)-1) % sampling_stride_x != 0:
    result += 1
  return result

def n_sampled_points_in_subdomain_y(subdomain_coordinate_y):
  result = (int)(np.ceil(n_fibers_in_subdomain_y(subdomain_coordinate_y) / sampling_stride_y))
  if subdomain_coordinate_y == n_subdomains_y-1 and (n_fibers_in_subdomain_y(subdomain_coordinate_y)-1) % sampling_stride_y != 0:
    result += 1
  return result

def n_sampled_points_in_subdomain_z(subdomain_coordinate_z):
  result = (int)(np.ceil(n_points_in_subdomain_z(subdomain_coordinate_z) / sampling_stride_z))
  if subdomain_coordinate_z == n_subdomains_z-1 and (n_points_in_subdomain_z(subdomain_coordinate_z)-1) % sampling_stride_z != 0:
    result += 1
  return result

# define 3D mesh 
# loop over nodes
node_positions = []
subdomain_coordinate_x = rank_no % n_subdomains_x
subdomain_coordinate_y = (int)(rank_no / n_subdomains_x) % n_subdomains_y
subdomain_coordinate_z = (int)(rank_no / n_subdomains_xy)

z_point_index_start = subdomain_coordinate_z * n_points_per_subdomain_z
z_point_index_end = z_point_index_start + n_points_in_subdomain_z(subdomain_coordinate_z)

# local number of elements
n_elements = [
    n_sampled_points_in_subdomain_x(subdomain_coordinate_x),
    n_sampled_points_in_subdomain_y(subdomain_coordinate_y), 
    n_sampled_points_in_subdomain_z(subdomain_coordinate_z)
  ]

n_global_points_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(n_subdomains_x)])
n_global_points_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(n_subdomains_y)])
n_global_points_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(n_subdomains_z)])

print("n_global_points: {} x {} x {} = sum({}) x sum({}) x sum({})".format(n_global_points_x, n_global_points_y, n_global_points_z, \
  [n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(n_subdomains_x)],\
  [n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(n_subdomains_y)],\
  [n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(n_subdomains_z)]))

# border subdomains have one element less than fibers
if subdomain_coordinate_x == n_subdomains_x-1:
  n_elements[0] -= 1
if subdomain_coordinate_y == n_subdomains_y-1:
  n_elements[1] -= 1
if subdomain_coordinate_z == n_subdomains_z-1:
  n_elements[2] -= 1

# loop over z point indices
for k in range(n_sampled_points_in_subdomain_z(subdomain_coordinate_z)):
  z_point_index = z_point_index_start + k*sampling_stride_z
  
  if subdomain_coordinate_z == n_subdomains_z-1 and k == n_sampled_points_in_subdomain_z(subdomain_coordinate_z)-1:
    z_point_index = z_point_index_end-1
    
  #print("{}: sampling_stride_z: {}, k: {}, z: {}/{}".format(rank_no, sampling_stride_z, k, z_point_index, z_point_index_end))
  
  # loop over fibers for own rank
  for j in range(n_sampled_points_in_subdomain_y(subdomain_coordinate_y)):
    fiber_in_subdomain_coordinate_y = j*sampling_stride_y
    
    if subdomain_coordinate_y == n_subdomains_y-1 and j == n_sampled_points_in_subdomain_y(subdomain_coordinate_y)-1:
      fiber_in_subdomain_coordinate_y = n_fibers_in_subdomain_y(subdomain_coordinate_y)-1
    
    #print("{}: sampling_stride_y: {}, j: {}, y: {}/{}".format(rank_no, sampling_stride_y, j, fiber_in_subdomain_coordinate_y, n_fibers_in_subdomain_y(subdomain_coordinate_y)))
      
    for i in range(n_sampled_points_in_subdomain_x(subdomain_coordinate_x)):
      fiber_in_subdomain_coordinate_x = i*sampling_stride_x
      
      if subdomain_coordinate_x == n_subdomains_x-1 and i == n_sampled_points_in_subdomain_x(subdomain_coordinate_x)-1:
        fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(subdomain_coordinate_x)-1
      
      #print("{}: sampling_stride_x: {}, i: {}, x: {}/{}".format(rank_no, sampling_stride_x, i, fiber_in_subdomain_coordinate_x, n_fibers_in_subdomain_x(subdomain_coordinate_x)))
      
      fiber_index = fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
        
      
      point = fibers[fiber_index][z_point_index]
      node_positions.append(point)
      #print("{}: {}".format(rank_no, point))
      
# on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
if subdomain_coordinate_x == n_subdomains_x-1:
  fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(subdomain_coordinate_x)-1


print("rank {}, n_elements: {}, subdomain coordinate ({},{},{})/({},{},{})".format(rank_no, n_elements, subdomain_coordinate_x, subdomain_coordinate_y, subdomain_coordinate_z, n_subdomains_x, n_subdomains_y, n_subdomains_z))
print("   fibers x: [{}, {}]".format(0, n_fibers_in_subdomain_x(subdomain_coordinate_x)))
print("   fibers y: [{}, {}]".format(0, n_fibers_in_subdomain_y(subdomain_coordinate_y)))
print("   points z: [{}, {}] ({})".format(z_point_index_start, z_point_index_end, n_points_in_subdomain_z(subdomain_coordinate_z)))
      
meshes["3Dmesh"] = {
  "nElements": n_elements,
  "nRanks": [n_subdomains_x, n_subdomains_y, n_subdomains_z],
  "nodePositions": node_positions,
  "inputMeshIsGlobal": False,
  "setHermiteDerivatives": False,
  "logKey": "3Dmesh"
}

# set Dirichlet BC values for bottom nodes to 0 and for top nodes to 1
potential_flow_dirichlet_bc = {}
for i in range(n_global_points_x*n_global_points_y):
  potential_flow_dirichlet_bc[i] = 0.0
  potential_flow_dirichlet_bc[(n_global_points_z-1)*n_global_points_x*n_global_points_y + i] = 1.0
    
print("potential_flow_dirichlet_bc: {}".format(potential_flow_dirichlet_bc))
    
if rank_no == 0:
  print("rank configuration: ")
  
  for subdomain_coordinate_y in range(n_subdomains_y):
    for subdomain_coordinate_x in range(n_subdomains_x):
      print("subdomain ({},{})".format(subdomain_coordinate_x, subdomain_coordinate_y))
      for rankNo in range(subdomain_coordinate_y*n_subdomains_x + subdomain_coordinate_x, n_ranks, n_subdomains_x*n_subdomains_y):
        print("rank {}".format(rankNo))
      print("  n_subdomains_z: {}".format(n_subdomains_z))
      for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
          print("  fiber {} in subdomain ({},{}) uses ranks {}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y), fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y, list(range(subdomain_coordinate_y*n_subdomains_x + subdomain_coordinate_x, n_ranks, n_subdomains_x*n_subdomains_y))))

config = {
  "scenarioName": scenario_name,
  "Meshes": meshes,
  "MappingsBetweenMeshes": {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(n_fibers_total)},
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": solver_type,
      "preconditionerType": "none",
    },
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
  "Coupling": {
    "timeStepWidth": dt_3D,  # 1e-1
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval" : 1,
    "endTime": end_time,
    "Term1": {      # monodomain fibers
      "MultipleInstances": {
        "nInstances": n_subdomains_xy,
        "instances": 
        [{
          "ranks": list(range(subdomain_coordinate_y*n_subdomains_x + subdomain_coordinate_x, n_ranks, n_subdomains_x*n_subdomains_y)),
          "StrangSplitting": {
            #"numberTimeSteps": 1,
            "timeStepWidth": dt_3D,  # 1e-1
            "logTimeStepWidthAsKey": "dt_3D",
            "durationLogKey": "duration_monodomain",
            "timeStepOutputInterval" : 1000,
            "endTime": dt_3D,

            "Term1": {      # CellML
              "MultipleInstances": {
                "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                "instances": 
                [{
                  "ranks": list(range(n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
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
                      "compilerFlags": "-fPIC -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared ",
                      #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
                      #"libraryFilename": "cellml_simd_lib.so",   # compiled library
                      "useGivenLibrary": False,
                      #"statesInitialValues": [],
                      #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
                      #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                      "setSpecificStatesCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "additionalArgument": i,
                      
                      "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
                      "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                      "meshName": "MeshFiber_{}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)),
                      "prefactor": 1.0,
                    },
                  },
                } for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y))],
              }
            },
            "Term2": {     # Diffusion
              "MultipleInstances": {
                "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                "instances": 
                [{
                  "ranks": list(range(n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "ImplicitEuler" : {
                    "initialValues": [],
                    #"numberTimeSteps": 1,
                    "timeStepWidth": dt_1D,  # 1e-5
                    "logTimeStepWidthAsKey": "dt_1D",
                    "durationLogKey": "duration_1D",
                    "timeStepOutputInterval": 1e4,
                    "dirichletBoundaryConditions": {0: -75.0036, -1: -75.0036},
                    "inputMeshIsGlobal": True,
                    "solverName": "implicitSolver",
                    "FiniteElementMethod" : {
                      "maxIterations": 1e4,
                      "relativeTolerance": 1e-10,
                      "inputMeshIsGlobal": True,
                      "meshName": "MeshFiber_{}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)),
                      "prefactor": Conductivity/(Am*Cm),
                      "solverName": "implicitSolver",
                    },
                    "OutputWriter" : [
                      #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fibre_"+str(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)), "binary": True, "fixedFormat": False, "combineFiles": True},
                      #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fibre_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                      #{"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
                      #{"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "binary":True, "onlyNodalValues":True},
                    ]
                  },
                } for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y))],
                "OutputWriter" : [
                  {"format": "Paraview", "outputInterval": int(1./dt_3D*output_timestep), "filename": "out/all_fibres", "binary": True, "fixedFormat": False, "combineFiles": True},
                ],
              },
            },
          }
        } for subdomain_coordinate_x in range(n_subdomains_x) for subdomain_coordinate_y in range(n_subdomains_y)]
      }
    },
    "Term2": {        # Bidomain
      "StaticBidomainSolver": {
        "timeStepWidth": dt_bidomain,
        "timeStepOutputInterval": 50,
        "solverName": "activationSolver",
        "inputIsGlobal": True,
        "PotentialFlow": {
          "FiniteElementMethod" : {  
            "meshName": "3Dmesh",
            "solverName": "potentialFlowSolver",
            "prefactor": 1.0,
            "dirichletBoundaryConditions": potential_flow_dirichlet_bc,
            "inputMeshIsGlobal": True,
          },
        },
        "Activation": {
          "FiniteElementMethod" : {  
            "meshName": "3Dmesh",
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
            ],
          },
        },
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./dt_3D*output_timestep), "filename": "out/3d", "binary": True, "fixedFormat": False, "combineFiles": False},
        ],
      }
    }
  }
}
