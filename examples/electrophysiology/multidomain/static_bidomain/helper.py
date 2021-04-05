# Multiple 1D fibers (monodomain) with 3D EMG (static bidomain), biceps geometry
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py
#
# if variables.fiber_file=cuboid.bin, it uses a small cuboid test example

import numpy as np
import pickle
import sys,os
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
#if rank_no == 0:
#  print("diffusion solver type: {}".format(variables.diffusion_solver_type))

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

if variables.n_subdomains != n_ranks:
  print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
  quit()
  

#############################

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

#############################
# create fat layer mesh
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

# parse parameters like number of points from the fat layer mesh file
parameters = []
for i in range(int(header_length/4.) - 1):
  int_raw = fat_mesh_file_handle.read(4)
  value = struct.unpack('i', int_raw)[0]
  parameters.append(value)

n_points_xy = parameters[0]
n_points_z = parameters[1]
n_points_x = (int)(np.sqrt(parameters[0]))
n_points_y = n_points_x

if "version 2" in header_str.decode("utf-8"):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
  n_points_x = parameters[2]
  n_points_y = parameters[3]

n_points = n_points_x*n_points_y*n_points_z

# output of global size
#if rank_no == 0:
#  print("    fat mesh, n points total:    {} ({} x {} x {})".format(n_points, n_points_x, n_points_y, n_points_z))
  
if False:
  # parse whole file, this is no longer needed as only the local node positions will get parsed
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

# determine number of nodes of the 3Dmesh
n_elements_3D_mesh_linear = variables.meshes["3Dmesh"]["nElements"]
n_points_3D_x = n_elements_3D_mesh_linear[0]
n_points_3D_y = n_elements_3D_mesh_linear[1]
n_points_3D_z = n_elements_3D_mesh_linear[2]

# if the own subdomain is at the (x+) boundary
if variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
  n_points_3D_x += 1

# if the own subdomain is at the (y+) boundary
if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
  n_points_3D_y += 1
  
# if the own subdomain is at the (z+) boundary
if variables.own_subdomain_coordinate_z == variables.n_subdomains_z - 1:
  n_points_3D_z += 1

if False:
  print("{}: 3Dmesh has {} local node positions ({} x {} x {} = {}), nElements: {}".format(rank_no, len(variables.meshes["3Dmesh"]["nodePositions"]), n_points_3D_x, n_points_3D_y, n_points_3D_z, n_points_3D_x*n_points_3D_y*n_points_3D_z, n_elements_3D_mesh_linear))

# determine nodes of the fat layer mesh
fat_mesh_node_indices = []

# range of points in z direction is [variables.z_point_index_start, variables.z_point_index_end)
# these variables were set by create_partitioned_meshes_for_settings
n_sampled_points_3D_in_own_subdomain_x = n_sampled_points_in_subdomain_x(variables.own_subdomain_coordinate_x)
n_sampled_points_3D_in_own_subdomain_y = n_sampled_points_in_subdomain_y(variables.own_subdomain_coordinate_y)
n_sampled_points_3D_in_own_subdomain_z = n_sampled_points_in_subdomain_z(variables.own_subdomain_coordinate_z)

n_points_on_previous_ranks_all_x = 0
n_points_on_previous_ranks_sampled_x = 0

if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
  #print("{}: previous subdomains: y+: {}".format(rank_no, list(range(variables.own_subdomain_coordinate_x))))
  n_points_on_previous_ranks_all_x = n_fibers_in_previous_subdomains_x(variables.own_subdomain_coordinate_x)
  n_points_on_previous_ranks_sampled_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.own_subdomain_coordinate_x)])
  

elif variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
  #print("{}: previous subdomains y+: {}, x+: {}".format(rank_no, list(range(variables.n_subdomains_x)), list(range(variables.n_subdomains_y-1, variables.own_subdomain_coordinate_y, -1))))
  n_points_on_previous_ranks_all_x = variables.n_fibers_x \
    + sum([n_fibers_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y-1, variables.own_subdomain_coordinate_y, -1)]) \
    - 1
  n_points_on_previous_ranks_sampled_x = variables.n_points_3D_mesh_global_x \
    + sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y-1, variables.own_subdomain_coordinate_y, -1)]) \
    - 1

n_points_on_previous_ranks_sampled_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.own_subdomain_coordinate_z)])
  
#print("{}: n_points_on_previous_ranks_all_x: {}".format(rank_no, n_points_on_previous_ranks_all_x))
#print("{}: z range: [{},{}), stride: {}".format(rank_no, variables.z_point_index_start, variables.z_point_index_end, variables.sampling_stride_z))

# loop over z point indices of the 3D mesh
for k in range(n_sampled_points_3D_in_own_subdomain_z):
  z_point_index = variables.z_point_index_start + k*variables.sampling_stride_z
  
  if variables.own_subdomain_coordinate_z == variables.n_subdomains_z-1 and k == n_sampled_points_3D_in_own_subdomain_z-1:
    z_point_index = variables.z_point_index_end-1
    
  #print("{}: sampling_stride_z: {}, k: {}, z: {}/{}".format(rank_no, variables.sampling_stride_z, k, z_point_index, variables.z_point_index_end))
  
  # loop over points of the fat layer mesh in y direction
  for j in range(0,n_points_y,variables.sampling_stride_fat):
    y_point_index = j
    
    # loop over the points in x direction of the fat layer mesh, this corresponds to x and negative y direction of the 3D mesh (see figure above in the code)
    
    # loop over points in 3D mesh in x direction
    x_point_index_offset = n_points_on_previous_ranks_all_x
    fat_mesh_n_points_x = 0
    if variables.own_subdomain_coordinate_y == variables.n_subdomains_y-1:
      
      
      for i_3D in range(n_sampled_points_3D_in_own_subdomain_x):
        x_point_index = n_points_on_previous_ranks_all_x + i_3D*variables.sampling_stride_x
        
        # on boundary rank set last node positions to be the boundary nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
        if variables.own_subdomain_coordinate_x == variables.n_subdomains_x-1 and i_3D == n_sampled_points_3D_in_own_subdomain_x-1:
          x_point_index = n_points_on_previous_ranks_all_x + n_fibers_in_subdomain_x(variables.own_subdomain_coordinate_x)-1
            
        # store index of node
        fat_mesh_node_indices.append([x_point_index, y_point_index, z_point_index])
        fat_mesh_n_points_x += 1
            
        #if j == 0 and k == 0:
        #  print("{}, at y+, x={}".format(rank_no, x_point_index, point))
            
      x_point_index_offset = x_point_index
      
    # loop over points in 3D mesh in negative y direction
    if variables.own_subdomain_coordinate_x == variables.n_subdomains_x-1:
        
      for j_3D in reversed(range(n_sampled_points_3D_in_own_subdomain_y)):
        y_index_3D_mesh = j_3D*variables.sampling_stride_y
        
        # on top boundary rank do not use the top node again, it was already visited within the x traversal
        if variables.own_subdomain_coordinate_y == variables.n_subdomains_y-1 and j_3D == n_sampled_points_3D_in_own_subdomain_y-1:
          continue
        
        x_point_index = x_point_index_offset + n_fibers_in_subdomain_y(variables.own_subdomain_coordinate_y)-1 - y_index_3D_mesh
        
        # store index of node
        fat_mesh_node_indices.append([x_point_index, y_point_index, z_point_index])
        fat_mesh_n_points_x += 1
        
        #if j == 0 and k == 0:
        #  print("{}, at x+, x={}".format(rank_no, x_point_index))

# load local nodes from file
fat_mesh_node_positions_local = []
for (index_x, index_y, index_z) in fat_mesh_node_indices:
  point = []
  # note, ordering in bin file is fastest in z, then in x then in y direction
  global_index = index_y*n_points_x*n_points_z + index_x*n_points_z + index_z
  
  # set file pointer to position of current index
  fat_mesh_file_handle.seek(32+10*4+global_index*3*8)
  
  # parse point
  for component_no in range(3):
    double_raw = fat_mesh_file_handle.read(8)
    value = struct.unpack('d', double_raw)[0]
    point.append(value)
  
  # store point in the list of local node positions of the fat mesh
  fat_mesh_node_positions_local.append(point)
  
# local size
n_sampled_points_y = len(range(0,n_points_y,variables.sampling_stride_fat))
fat_mesh_n_points = [fat_mesh_n_points_x, n_sampled_points_y, n_sampled_points_3D_in_own_subdomain_z]
fat_mesh_n_elements = [fat_mesh_n_points[0]-1, fat_mesh_n_points[1]-1, fat_mesh_n_points[2]-1]

# regarding x direction, if in interior of fat mesh (i.e., not bottom right subdomain), adjust number of elements in x direction (which is in negative y direction)
if not (variables.own_subdomain_coordinate_y == 0 and variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1):
  fat_mesh_n_elements[0] += 1
# regarding z direction, if in interior of fat mesh, adjust number of elements in z direction
if variables.own_subdomain_coordinate_z != variables.n_subdomains_z - 1:
  fat_mesh_n_elements[2] += 1
          
# store values to be used in postprocess callback function
variables.fat_mesh_n_points_local = fat_mesh_n_points
variables.fat_mesh_n_points_global = [variables.n_points_3D_mesh_global_x+variables.n_points_3D_mesh_global_y-1, n_sampled_points_y, variables.n_points_3D_mesh_global_z]
variables.fat_mesh_index_offset = [n_points_on_previous_ranks_sampled_x, 0, n_points_on_previous_ranks_sampled_z]

# debugging output
if False:
  # if the own subdomain is at the (y+) boundary
  if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
    # output top row node positions of 3D mesh
    for i in range(n_points_3D_x):
      point = variables.meshes["3Dmesh"]["nodePositions"][0*n_points_3D_x*n_points_3D_y + (n_points_3D_y-1)*n_points_3D_x + i]
      
      print("{}: 3Dmesh top node i={}, {}".format(rank_no, i, point))
    
  # if the own subdomain is at the (x+) boundary
  if variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
    # output top row node positions of 3D mesh
    for j in reversed(range(n_points_3D_y)):
      point = variables.meshes["3Dmesh"]["nodePositions"][0*n_points_3D_x*n_points_3D_y + j*n_points_3D_x + (n_points_3D_x-1)]
      
      print("{}: 3Dmesh right node j={}, {}".format(rank_no, j, point))
    
  # output fat layer mesh on bottom
  for i in range(fat_mesh_n_points[0]):
    point = fat_mesh_node_positions_local[i]
    print("{}: fat mesh bottom node i={}, {}".format(rank_no, i, point))
  
fat_mesh_n_ranks = [variables.n_subdomains_x + variables.n_subdomains_y - 1, 1, variables.n_subdomains_z]

if False:
  print("{}: Rank {} has subset of fat layer mesh: {} x {} x {} = {} = {} nodes, n elements local: {}, n points local: {}, global: {}".format(rank_no, rank_no, \
    fat_mesh_n_points[0], fat_mesh_n_points[1], fat_mesh_n_points[2], fat_mesh_n_points[0]*fat_mesh_n_points[1]*fat_mesh_n_points[2], 
    len(fat_mesh_node_positions_local), fat_mesh_n_elements, variables.fat_mesh_n_points_local, variables.fat_mesh_n_points_global))

# determine all ranks that participate in computing the mesh (global nos)
variables.fat_global_rank_nos = []
for coordinate_z in range(variables.n_subdomains_z):
  for coordinate_x in range(variables.n_subdomains_x):
    fat_rank_no = coordinate_z * variables.n_subdomains_xy + (variables.n_subdomains_y-1)*variables.n_subdomains_x + coordinate_x
    variables.fat_global_rank_nos.append(fat_rank_no)
    
  for coordinate_y in range(variables.n_subdomains_y-2,-1,-1):
    fat_rank_no = coordinate_z * variables.n_subdomains_xy + coordinate_y*variables.n_subdomains_x + (variables.n_subdomains_x-1)
    variables.fat_global_rank_nos.append(fat_rank_no)
    
if False:
  print("{}: Fat mesh will be computed by the following ranks: fat_global_rank_nos: {}, fat_mesh_n_ranks: {}".
    format(rank_no, variables.fat_global_rank_nos, fat_mesh_n_ranks))

if rank_no == 0:
  print("    fat mesh, n points total:    {} ({} x {} x {}), (per process: {} x {} x {} = {})".format(variables.fat_mesh_n_points_global[0]*variables.fat_mesh_n_points_global[1]*variables.fat_mesh_n_points_global[2], variables.fat_mesh_n_points_global[0], variables.fat_mesh_n_points_global[1], variables.fat_mesh_n_points_global[2], fat_mesh_n_points[0], fat_mesh_n_points[1], fat_mesh_n_points[2], fat_mesh_n_points[0]*fat_mesh_n_points[1]*fat_mesh_n_points[2]))
  
#print("fat mesh:")
#print(fat_mesh_node_positions_local)

variables.meshes["3DFatMesh"] = {
  "nElements": fat_mesh_n_elements,
  "nRanks": fat_mesh_n_ranks,
  "rankNos": variables.fat_global_rank_nos,
  "nodePositions": fat_mesh_node_positions_local,
  "inputMeshIsGlobal": False,
  "setHermiteDerivatives": False,
  "logKey": "3DFatMesh"
}

# create mappings between meshes
variables.mappings_between_meshes = {
  "3Dmesh": {
    "name": "3DFatMesh", 
    "xiTolerance": 1e-2,
    "enableWarnings": True,
    "compositeUseOnlyInitializedMappings": False,
    "fixUnmappedDofs": False,
    "defaultValue": 0,
  }
}
# only include overlapping elements

# set output writer    
variables.output_writer_fibers = []
variables.output_writer_elasticity = []
variables.output_writer_emg = []
variables.output_writer_0D_states = []

subfolder = ""
if variables.paraview_output:
  if variables.adios_output:
    subfolder = "paraview/"
  variables.output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  if variables.states_output:
    variables.output_writer_0D_states.append({"format": "Paraview", "outputInterval": 1, "filename": "out/" + subfolder + variables.scenario_name + "/0D_states", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})

if variables.adios_output:
  if variables.paraview_output:
    subfolder = "adios/"
  variables.output_writer_emg.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "useFrontBackBuffer": False, "fileNumbering": "incremental"})
  #variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/fibers", "combineNInstances": 1, "useFrontBackBuffer": False}

if variables.python_output:
  if variables.adios_output:
    subfolder = "python/"
  variables.output_writer_emg.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fileNumbering": "incremental"})
  
if variables.exfile_output:
  if variables.adios_output:
    subfolder = "exfile/"
  variables.output_writer_emg.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "fileNumbering": "incremental"})
  variables.output_writer_elasticity.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "fileNumbering": "incremental"})
  
# set variable mappings for cellml model
if "hodgkin_huxley" in variables.cellml_file:
  # parameters: I_stim
  variables.mappings = {
    ("parameter", 0):           ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("connectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0]                         # initial value for stimulation current
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

elif "shorten" in variables.cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  variables.mappings = {
    ("parameter", 0):           ("algebraic", "wal_environment/I_HH"), # parameter is algebraic 32
    ("parameter", 1):           ("constant", "razumova/L_x"),             # parameter is constant 65, fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension
    ("connectorSlot", 0): ("state", "wal_environment/vS"),          # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                        # stimulation current I_stim, fiber stretch λ
  variables.nodal_stimulation_current = 1200.                             # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "slow_TK_2014" in variables.cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
  # parameters: I_stim, fiber stretch λ
  variables.mappings = {
    ("parameter", 0):           ("constant", "wal_environment/I_HH"), # parameter 0 is constant 54 = I_stim
    ("parameter", 1):           ("constant", "razumova/L_S"),         # parameter 1 is constant 67 = fiber stretch λ
    ("connectorSlot", 0): ("state", "wal_environment/vS"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "razumova/stress"),  # expose algebraic 12 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                    # wal_environment/I_HH = I_stim, razumova/L_S = λ
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_2016_08_22" in variables.cellml_file :   # this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 8 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/velo"),        # parameter 2 is constant 9 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/sigma"),   # expose algebraic 0 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_Titin" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 11 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/rel_velo"),    # parameter 2 is constant 12 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/ActiveStress"),   # expose algebraic 4 = γ to the operator splitting
    ("connectorSlot", 2): ("algebraic", "Razumova/Activation"),     # expose algebraic 5 = α to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

# --------------------------
# callback functions

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
    
