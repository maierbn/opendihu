# Multiple 1D fibers (monodomain) with 3D EMG (static bidomain), biceps geometry
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py
#
# if variables.fiber_file=cuboid.bin, it uses a small cuboid test example

import numpy as np
import pickle
import sys
import struct
import argparse
sys.path.insert(0, '..')
import variables    # file variables.py
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# generate cuboid fiber file
if "cuboid.bin" in variables.fiber_file:
  
  if variables.n_fibers_y is None:
    variables.n_fibers_x = 4
    variables.n_fibers_y = variables.n_fibers_x
    variables.n_points_whole_fiber = 20
  
  size_x = variables.n_fibers_x * 0.1
  size_y = variables.n_fibers_y * 0.1
  size_z = variables.n_points_whole_fiber / 100.
  
  if rank_no == 0:
    print("create cuboid.bin with size [{},{},{}], n points [{},{},{}]".format(size_x, size_y, size_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber))
    
    # write header
    with open(variables.fiber_file, "wb") as outfile:
      
      # write header
      header_str = "opendihu self-generated cuboid  "
      outfile.write(struct.pack('32s',bytes(header_str, 'utf-8')))   # 32 bytes
      outfile.write(struct.pack('i', 40))  # header length
      outfile.write(struct.pack('i', variables.n_fibers_x*variables.n_fibers_y))   # n_fibers
      outfile.write(struct.pack('i', variables.n_points_whole_fiber))   # variables.n_points_whole_fiber
      outfile.write(struct.pack('i', 0))   # nBoundaryPointsXNew
      outfile.write(struct.pack('i', 0))   # nBoundaryPointsZNew
      outfile.write(struct.pack('i', 0))   # nFineGridFibers_
      outfile.write(struct.pack('i', 1))   # nRanks
      outfile.write(struct.pack('i', 1))   # nRanksZ
      outfile.write(struct.pack('i', 0))   # nFibersPerRank
      outfile.write(struct.pack('i', 0))   # date
    
      # loop over points
      for y in range(variables.n_fibers_y):
        for x in range(variables.n_fibers_x):
          for z in range(variables.n_points_whole_fiber):
            point = [x*(float)(size_x)/(variables.n_fibers_x), y*(float)(size_y)/(variables.n_fibers_y), z*(float)(size_z)/(variables.n_points_whole_fiber)]
            outfile.write(struct.pack('3d', point[0], point[1], point[2]))   # data point

# output diffusion solver type
if rank_no == 0:
  print("diffusion solver type: {}".format(variables.diffusion_solver_type))

variables.load_fiber_data = True   # load all local node positions from fiber_file, in order to infer partitioning for fat_layer mesh

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    variables.fiber_file, variables.load_fiber_data,
    variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z, variables.generate_linear_3d_mesh, variables.generate_quadratic_3d_mesh)
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result
  
variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# set fat layer mesh
# load file, the file is assumed to be small enough to be loaded completely by all ranks

try:
  fat_mesh_file_handle = open(variables.fat_mesh_file, "rb")
except:
  print("Error: Could not open fat mesh file \"{}\"".format(variables.fat_mesh_file))
  quit()

# parse file header to extract mesh dimensions
bytes_raw = fat_mesh_file_handle.read(32)
header_str = struct.unpack('32s', bytes_raw)[0]
header_length_raw = fat_mesh_file_handle.read(4)
header_length = struct.unpack('i', header_length_raw)[0]

# parse parameters in the file
parameters = []
for i in range(int(header_length/4.) - 1):
  double_raw = fat_mesh_file_handle.read(4)
  value = struct.unpack('i', double_raw)[0]
  parameters.append(value)

n_points_xy = parameters[0]
n_points_z = parameters[1]
n_points_x = (int)(np.sqrt(parameters[0]))
n_points_y = n_points_x

if "version 2" in header_str.decode("utf-8"):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
  n_points_x = parameters[2]
  n_points_y = parameters[3]

n_points = n_points_x*n_points_y*n_points_z

if rank_no == 0:
  print("fat mesh, n points:    {} ({} x {} x {})".format(n_points, n_points_x, n_points_y, n_points_z))
  
# parse whole file
fat_mesh_node_positions = [None for _ in range(n_points)]
for j in range(n_points_y):
  for i in range(n_points_x):
    for k in range(n_points_z):
      
      point = []
      for component_no in range(3):
        double_raw = fat_mesh_file_handle.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
      fat_mesh_node_positions[k*n_points_xy + j*n_points_x + i] = point

# The fat mesh "3DFatMesh" touches the 3Dmesh of intramuscular EMG at an interface surface.
# Do the domain decomposition such that neighbouring elements across this interface are on the same rank.
# This means that the fat mesh will not be computed by all ranks.

#     /xxxxxx
# ^  /xxxx/xx     x  = 3D fat layer mesh
# |  |_|_|/xx    |_| = intramuscular 3D mesh
# y  |_|_|/x/
#  x-->
#

fat_mesh_node_positions_local = []

n_elements_3D_mesh_linear = variables.meshes["3Dmesh"]["nElements"]
n_points_3D_x = (n_elements_3D_mesh_linear[0]+1)
n_points_3D_y = (n_elements_3D_mesh_linear[1]+1)
n_points_3D_z = (n_elements_3D_mesh_linear[2]+1)

print("3Dmesh has {} node positions ({} x {} x {} = {})".format(len(variables.meshes["3Dmesh"]["nodePositions"]), n_points_3D_x, n_points_3D_y, n_points_3D_z, n_points_3D_x*n_points_3D_y*n_points_3D_z))

pz0 = variables.meshes["3Dmesh"]["nodePositions"][0*n_points_3D_x*n_points_3D_y + 0*n_points_3D_x + 0]
pz1 = variables.meshes["3Dmesh"]["nodePositions"][(n_points_3D_z-1)*n_points_3D_x*n_points_3D_y + 0*n_points_3D_x + 0]

# find indices k in z direction of the local part of the fat layer mesh
index_k_start = None
index_k_end = None

min_distance_start = None
min_distance_end = None

# loop over z indices of fat_mesh_node_positions
for k in range(n_points_z):
  p = fat_mesh_node_positions[k*n_points_xy + 0*n_points_x + 0]
  
  distance_pz0 = np.linalg.norm(np.array(p) - np.array(pz0)) 
  distance_pz1 = np.linalg.norm(np.array(p) - np.array(pz1)) 
  
  if min_distance_start is None or distance_pz0 < min_distance_start:
    min_distance_start = distance_pz0
    index_k_start = k
    
  if min_distance_end is None or distance_pz1 < min_distance_end:
    min_distance_end = distance_pz1
    index_k_end = k
  
index_i_start = None
index_i_end = None
    
# if the own subdomain is at the (y+) boundary (top in sketch)
if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
  
  px0 = variables.meshes["3Dmesh"]["nodePositions"][(n_points_3D_y-1)*n_points_3D_x + 0]
  px1 = variables.meshes["3Dmesh"]["nodePositions"][(n_points_3D_y-1)*n_points_3D_x + n_points_3D_x-1]
  
  # find range [index_i_start, index_i_end] for x of local node positions of fat mesh
  index_i_start = None
  index_i_end = None
  min_distance_start = None
  min_distance_end = None
  
  for i in range(n_points_x):
    p = fat_mesh_node_positions[0*n_points_x + i]
    
    distance_px0 = np.linalg.norm(np.array(p) - np.array(px0)) 
    distance_px1 = np.linalg.norm(np.array(p) - np.array(px1)) 
    
    if min_distance_start is None or distance_px0 < min_distance_start:
      min_distance_start = distance_px0
      index_i_start = i
      
    if min_distance_end is None or distance_px1 < min_distance_end:
      min_distance_end = distance_px1
      index_i_end = i

# if the own subdomain is at the (x+) boundary (right in sketch)
if variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
  
  py0 = variables.meshes["3Dmesh"]["nodePositions"][(n_points_3D_y-1)*n_points_3D_x + n_points_3D_x-1]
  py1 = variables.meshes["3Dmesh"]["nodePositions"][0*n_points_3D_x + n_points_3D_x-1]
  
  # find range [index_i_start, index_i_end] for x of local node positions of fat mesh
  min_distance_start = None
  min_distance_end = None
  
  for i in range(n_points_x):
    p = fat_mesh_node_positions[0*n_points_x + i]
    
    distance_py0 = np.linalg.norm(np.array(p) - np.array(py0)) 
    distance_py1 = np.linalg.norm(np.array(p) - np.array(py1)) 
    
    if not variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
      if min_distance_start is None or distance_py0 < min_distance_start:
        min_distance_start = distance_py0
        index_i_start = i
      
    if min_distance_end is None or distance_py1 < min_distance_end:
      min_distance_end = distance_py1
      index_i_end = i

# add all local nodes
previous_point = None
previous_info = None
for k in range(index_k_start,index_k_end+1):
  
  # do not include ghost nodes (if not at z+ boundary)
  if k == index_k_end and index_k_end != n_points_z-1:
    continue
  
  for j in range(n_points_y):
    for i in range(index_i_start,index_i_end+1):
          
      # do not include ghost nodes (if not at x+ boundary)
      if i == index_i_end and index_i_end != n_points_x-1:
        continue
        
      point = fat_mesh_node_positions[k*n_points_xy + j*n_points_x + i]
      fat_mesh_node_positions_local.append(point)
      
      if previous_point is not None and np.linalg.norm(np.array(previous_point) - np.array(point)) < 1e-8:
        print("Error, twice the same point! At i,j,k={},{},{} point {} (previous: {})".format(i,j,k,point, previous_info))
      
      previous_point = point
      previous_info = [i,j,k]

# local size
fat_mesh_n_points = [index_i_end+1 - index_i_start, n_points_y, index_k_end+1 - index_k_start]
fat_mesh_n_elements = [fat_mesh_n_points[0]-1, fat_mesh_n_points[1]-1, fat_mesh_n_points[2]-1]

# regarding x direction, if in interior
if index_i_end != n_points_x-1:
  fat_mesh_n_points[0] -= 1
# regarding z direction, if in interior
if index_k_end != n_points_z-1:
  fat_mesh_n_points[2] -= 1

fat_mesh_n_ranks = [variables.n_subdomains_x + variables.n_subdomains_y - 1, 1, variables.n_subdomains_z]

variables.fat_dirichlet_bc = {}
for k in range(fat_mesh_n_points[2]):
  for i in range(fat_mesh_n_points[0]):
    j = 0
    dof_no_local = k * fat_mesh_n_points[0]*fat_mesh_n_points[1] + j * fat_mesh_n_points[0] + i
    variables.fat_dirichlet_bc[dof_no_local] = dof_no_local

print("Fat mesh on rank {}, subset i: [{},{}], k: [{},{}], {} x {} x {} = {} = {} nodes".format(rank_no, index_i_start, index_i_end, index_k_start, index_k_end, \
  fat_mesh_n_points[0], fat_mesh_n_points[1], fat_mesh_n_points[2], fat_mesh_n_points[0]*fat_mesh_n_points[1]*fat_mesh_n_points[2], len(fat_mesh_node_positions_local) ))
#print("dofs: ",variables.fat_dirichlet_bc)


# determine all ranks that participate in computing the mesh (global nos)
variables.fat_global_rank_nos = []
for coordinate_z in range(variables.n_subdomains_z):
  for coordinate_x in range(variables.n_subdomains_x):
    rank_no = coordinate_z * variables.n_subdomains_xy + (variables.n_subdomains_y-1)*variables.n_subdomains_x + coordinate_x
    variables.fat_global_rank_nos.append(rank_no)
    
  for coordinate_y in range(variables.n_subdomains_y-2,-1,-1):
    rank_no = coordinate_z * variables.n_subdomains_xy + coordinate_y*variables.n_subdomains_x + (variables.n_subdomains_x-1)
    variables.fat_global_rank_nos.append(rank_no)


print("Fat mesh on rank {}, fat_global_rank_nos: {}, fat_mesh_n_ranks: {}".format(rank_no, variables.fat_global_rank_nos, fat_mesh_n_ranks))

#print("fat mesh:")
#print(fat_mesh_node_positions_local)

variables.meshes["3DFatMesh"] = {
  "nElements": fat_mesh_n_elements,
  "nRanks": fat_mesh_n_ranks,
  "nodePositions": fat_mesh_node_positions_local,
  "inputMeshIsGlobal": False,
  "setHermiteDerivatives": False,
  "logKey": "3DFatMesh"
}

if False:
  print("settings 3DFatMesh: ")
  with open("3DFatMesh","w") as f:
    f.write(str(variables.meshes["3DFatMesh"]))

# create mappings between meshes
variables.mappings_between_meshes = {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(variables.n_fibers_total)}
variables.mappings_between_meshes.update({"3Dmesh": {"name": "3DFatMesh", "xiTolerance": 1e-2, "defaultValue": 0}})    # only include overlapping elements

# set output writer    
variables.output_writer_fibers = []
variables.output_writer_elasticity = []
variables.output_writer_emg = []
variables.output_writer_0D_states = []

subfolder = ""
if variables.paraview_output:
  if variables.adios_output:
    subfolder = "paraview/"
  variables.output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fixedFormat": False, "combineFiles": True})
  variables.output_writer_elasticity.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fixedFormat": False, "combineFiles": True})
  variables.output_writer_fibers.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True, "fixedFormat": False, "combineFiles": True})
  if variables.states_output:
    variables.output_writer_0D_states.append({"format": "Paraview", "outputInterval": 1, "filename": "out/" + subfolder + variables.scenario_name + "/0D_states", "binary": True, "fixedFormat": False, "combineFiles": True})

if variables.adios_output:
  if variables.paraview_output:
    subfolder = "adios/"
  variables.output_writer_emg.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "useFrontBackBuffer": False, "combineNInstances": 1})
  variables.output_writer_elasticity.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "useFrontBackBuffer": False})
  variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "combineNInstances": variables.n_subdomains_xy, "useFrontBackBuffer": False})
  #variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/fibers", "combineNInstances": 1, "useFrontBackBuffer": False}

if variables.python_output:
  if variables.adios_output:
    subfolder = "python/"
  variables.output_writer_emg.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True})
  variables.output_writer_elasticity.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True})
  variables.output_writer_fibers.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True})

if variables.exfile_output:
  if variables.adios_output:
    subfolder = "exfile/"
  variables.output_writer_emg.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg"})
  variables.output_writer_elasticity.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity"})
  variables.output_writer_fibers.append({"format": "Exfile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/fibers"})

# set values for cellml model
if "shorten" in variables.cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  variables.parameters_used_as_algebraic = [32]    # 
  variables.parameters_used_as_constant = [65]        # fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension. CONSTANTS[65] in the shorten model
  variables.parameters_initial_values = [0.0, 1.0]    # stimulation current I_stim, fiber stretch λ, OpenCMISS generated files: OC_KNOWN will be set by this
  variables.nodal_stimulation_current = 1200.
  variables.output_state_index = 0                    # use state 0 = Vm
  variables.output_algebraic_index = []            # do not use any algebraic
  
elif "hodgkin_huxley" in variables.cellml_file:
  # parameters: I_stim
  variables.parameters_used_as_algebraic = []
  variables.parameters_used_as_constant = [2]
  variables.parameters_initial_values = [0.0]
  variables.nodal_stimulation_current = 40.
  variables.output_state_index = 0                    # use state 0 = Vm
  variables.output_algebraic_index = []            # do not use any algebraic

elif "slow_TK_2014" in variables.cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
  # parameters: I_stim, fiber stretch λ
  variables.parameters_used_as_algebraic = []
  variables.parameters_used_as_constant = [54, 67]     # wal_environment/I_HH = I_stim, razumova/L_S = λ
  variables.parameters_initial_values = [0.0, 1.0]     # wal_environment/I_HH = I_stim, razumova/L_S = λ
  variables.nodal_stimulation_current = 40. 
  variables.output_state_index = 0                     # use state 0, wal_environment/vS = Vm
  variables.output_algebraic_index = 12             # use algebraic 12, razumova/stress = γ
  
elif "Aliev_Panfilov_Razumova_2016_08_22" in variables.cellml_file :   # this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.parameters_used_as_algebraic = []
  variables.parameters_used_as_constant = [0, 8, 9]    # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
  variables.parameters_initial_values = [0, 1, 0]      # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
  variables.nodal_stimulation_current = 40. 
  variables.output_state_index = 0                     # use state 0, Aliev_Panfilov/V_m = Vm
  variables.output_algebraic_index = 0              # no algebraics are used
  
elif "Aliev_Panfilov_Razumova_Titin" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.parameters_used_as_algebraic = []
  variables.parameters_used_as_constant = [0, 11, 12]  # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
  variables.parameters_initial_values = [0, 1, 0]      # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
  variables.nodal_stimulation_current = 40. 
  variables.output_state_index = 0                     # use state 0, Aliev_Panfilov/V_m = Vm
  variables.output_algebraic_index = [4,5]          # Razumova/ActiveStress = γ, Razumova/Activation = α 
  

# callback functions
# --------------------------
def get_motor_unit_no(fiber_no):
  return int(variables.fiber_distribution[fiber_no % len(variables.fiber_distribution)]-1)

def get_diffusion_prefactor(fiber_no, mu_no):
  return variables.get_conductivity(fiber_no, mu_no) / (variables.get_am(fiber_no, mu_no) * variables.get_cm(fiber_no, mu_no))

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(fiber_no)*alpha)
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(variables.firing_times,0)
  
  #if variables.firing_times[index % n_firing_times, mu_no] == 1:
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), variables.firing_times[index % n_firing_times, mu_no], "true" if variables.firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return variables.firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fiber_no):
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = variables.innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)
  
  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = variables.nodal_stimulation_current
  else:
    stimulation_current = 0.
  
  first_dof_global = dof_nos_global[0]
  last_dof_global = dof_nos_global[-1]
    
  for node_no_global in nodes_to_stimulate_global:
    if first_dof_global <= node_no_global <= last_dof_global:
      # get local no for global no (1D)
      dof_no_local = node_no_global - first_dof_global
      parameters[dof_no_local] = stimulation_current
 
# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fiber_no):
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = variables.innervation_zone_width*100  # 100 nodes per cm
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

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time)) 
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:  
    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = variables.innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-1)
    if innervation_node_global < n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+1)
    if rank_no == 0:
      print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)

# load MU distribution and firing times
variables.fiber_distribution = np.genfromtxt(variables.fiber_distribution_file, delimiter=" ")
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
  print("Debugging output about fiber firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(variables.firing_times,0)
  for fiber_no_index in range(variables.n_fibers_total):
    if fiber_no_index % 100 == 0:
      t_algebraic = timeit.default_timer()
      if t_algebraic - t_start > 100:
        print("Note: break after {}/{} fibers ({:.0f}%) because it already took {:.3f}s".format(fiber_no_index,variables.n_fibers_total,100.0*fiber_no_index/(variables.n_fibers_total-1.),t_algebraic - t_start))
        break
    
    first_stimulation = None
    for current_time in np.linspace(0,1./variables.stimulation_frequency*n_firing_times,n_firing_times):
      if fiber_gets_stimulated(fiber_no_index, variables.stimulation_frequency, current_time):
        first_stimulation = current_time
        break
    mu_no = get_motor_unit_no(fiber_no_index)
    first_stimulation_info.append([fiber_no_index,mu_no,first_stimulation])
  
  first_stimulation_info.sort(key=lambda x: 1e6+1e-6*x[1]+1e-12*x[0] if x[2] is None else x[2]+1e-6*x[1]+1e-12*x[0])
  
  print("First stimulation times")
  print("    Time  MU fibers")
  n_stimulated_mus = 0
  n_not_stimulated_mus = 0
  stimulated_fibers = []
  last_time = 0
  last_mu_no = first_stimulation_info[0][1]
  for stimulation_info in first_stimulation_info:
    mu_no = stimulation_info[1]
    fiber_no = stimulation_info[0]
    if mu_no == last_mu_no:
      stimulated_fibers.append(fiber_no)
    else:
      if last_time is not None:
        if len(stimulated_fibers) > 10:
          print("{:8.2f} {:3} {} (only showing first 10, {} total)".format(last_time,last_mu_no,str(stimulated_fibers[0:10]),len(stimulated_fibers)))
        else:
          print("{:8.2f} {:3} {}".format(last_time,last_mu_no,str(stimulated_fibers)))
        n_stimulated_mus += 1
      else:
        if len(stimulated_fibers) > 10:
          print("  never stimulated: MU {:3}, fibers {} (only showing first 10, {} total)".format(last_mu_no,str(stimulated_fibers[0:10]),len(stimulated_fibers)))
        else:
          print("  never stimulated: MU {:3}, fibers {}".format(last_mu_no,str(stimulated_fibers)))
        n_not_stimulated_mus += 1
      stimulated_fibers = [fiber_no]

    last_time = stimulation_info[2]
    last_mu_no = mu_no
    
  print("stimulated MUs: {}, not stimulated MUs: {}".format(n_stimulated_mus,n_not_stimulated_mus))

  t_end = timeit.default_timer()
  print("duration of assembling this list: {:.3f} s\n".format(t_end-t_start))  
  
# compute partitioning
if rank_no == 0:
  if n_ranks != variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z:
    print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    quit()
  
variables.n_fibers_per_subdomain_x = (int)(variables.n_fibers_x / variables.n_subdomains_x)
variables.n_fibers_per_subdomain_y = (int)(variables.n_fibers_y / variables.n_subdomains_y)
variables.n_points_per_subdomain_z = (int)(variables.n_points_whole_fiber / variables.n_subdomains_z)

####################################
# set Dirichlet BC for the flow problem

n_points_3D_mesh_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_global = n_points_3D_mesh_global_x*n_points_3D_mesh_global_y*n_points_3D_mesh_global_z
 
# set Dirichlet BC values for bottom nodes to 0 and for top nodes to 1
variables.potential_flow_dirichlet_bc = {}
for i in range(n_points_3D_mesh_global_x*n_points_3D_mesh_global_y):
  variables.potential_flow_dirichlet_bc[i] = 0.0
  variables.potential_flow_dirichlet_bc[(n_points_3D_mesh_global_z-1)*n_points_3D_mesh_global_x*n_points_3D_mesh_global_y + i] = 1.0
    
# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top
variables.use_elasticity_dirichlet_bc = {}
for i in range(n_points_3D_mesh_global_x*n_points_3D_mesh_global_y):
  variables.use_elasticity_dirichlet_bc[(n_points_3D_mesh_global_z-1)*n_points_3D_mesh_global_x*n_points_3D_mesh_global_y + i] = 0.0
    
# Neumann BC at bottom nodes, traction downwards
nx = n_points_3D_mesh_global_x-1
ny = n_points_3D_mesh_global_y-1
nz = n_points_3D_mesh_global_z-1
variables.use_elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": [0.0,0.0,-1.0], "face": "2-"} for j in range(ny) for i in range(nx)]
#variables.use_elasticity_neumann_bc = []
    
# sanity checking at the end, is disabled and can be copied to after the config in the real settings file
if False:
  # check coupling instances
  multiple_instances = config["Coupling"]["Term1"]["MultipleInstances"]
  n_instances = multiple_instances["nInstances"]
  instances_size = len(multiple_instances["instances"])
  print("n_instances: {}".format(n_instances))

  print("n subdomains: {} x {}".format(variables.n_subdomains_x, variables.n_subdomains_y))
  print("n_fibers_per_subdomain_x: {} {}".format(variables.n_fibers_per_subdomain_x, variables.n_fibers_per_subdomain_y))

  for subdomain_coordinate_y in range(variables.n_subdomains_y):
    print("n_fibers_in_subdomain_y({}) = {}".format(subdomain_coordinate_y, n_fibers_in_subdomain_y(subdomain_coordinate_y)))

  print("--")
  for subdomain_coordinate_x in range(variables.n_subdomains_x):
    print("n_fibers_in_subdomain_x({}) = {}".format(subdomain_coordinate_x, n_fibers_in_subdomain_x(subdomain_coordinate_x)))
  print("--")

  # check fiber no
  counter = 0
  for subdomain_coordinate_y in range(variables.n_subdomains_y):
    for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
      for subdomain_coordinate_x in range(variables.n_subdomains_x):
        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
          no = get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
          if no != counter:
            print("error: get_fiber_no({},{},{},{}) = {}, counter = {}".format(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,no,counter))
          else:
            print("   ok: get_fiber_no({},{},{},{}) = {}, counter = {}".format(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,no,counter))
          counter += 1
          

  if n_instances != instances_size or instances_size == 0:
    print("Error with top-level multiple instances: nInstances: {}, size of instances: {}".format(n_instances, instances_size))

  # loop over inner instances
  for i in range(n_instances):
    multiple_instances0 = multiple_instances["instances"][i]["StrangSplitting"]["Term1"]["MultipleInstances"]
    n_instances0 = multiple_instances0["nInstances"]
    instances_size0 = len(multiple_instances0["instances"])
    
    if n_instances0 != instances_size0 or instances_size0 == 0:
      print("Error with Term1 {} multiple instances: nInstances: {}, size of instances: {}".format(i ,n_instances0, instances_size0))
    
    multiple_instances1 = multiple_instances["instances"][i]["StrangSplitting"]["Term2"]["MultipleInstances"]
    n_instances1 = multiple_instances1["nInstances"]
    instances_size1 = len(multiple_instances1["instances"])

    if n_instances1 != instances_size1 or instances_size1 == 0:
      print("Error with Term2 {} multiple instances: nInstances: {}, size of instances: {}".format(i, n_instances1, instances_size1))
