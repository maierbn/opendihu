# Multidomain helper script
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py

import numpy as np
import scipy
import pickle
import sys,os
import struct
import argparse
import time
sys.path.insert(0, '..')
import variables    # file variables.py
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

if variables.n_subdomains != n_ranks:
  print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\033[0m\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
  quit()

variables.relative_factors_file = "compartments_relative_factors.{}.{}_mus_stride_{}x{}x{}_partitioning_{}x{}x{}".format(
  os.path.basename(variables.fiber_file),
  len(variables.motor_units),
  variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z,
  variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z,
)

include_global_node_positions = False
if not os.path.exists(variables.relative_factors_file) and rank_no == 0:
  include_global_node_positions = True

  
#############################
# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    variables.fiber_file, variables.load_fiber_data,
    variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z, variables.generate_linear_3d_mesh, variables.generate_quadratic_3d_mesh,
    fiber_set_rank_nos=False, have_fibers=False, include_global_node_positions=include_global_node_positions)
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result
  
variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

variables.n_compartments = len(variables.motor_units)

#############################
# create fat layer mesh
try:
  fat_mesh_file_handle = open(variables.fat_mesh_file, "rb")
except:
  print("\033[0;31mError: Could not open fat mesh file \"{}\"\033[0m".format(variables.fat_mesh_file))
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

fat_mesh_quadratic_n_elements = [
  int(fat_mesh_n_elements[0]/2),
  int(fat_mesh_n_elements[1]/2),
  int(fat_mesh_n_elements[2]/2)
]

#print("fat_mesh_n_elements:           {}".format(fat_mesh_n_elements))
#print("fat_mesh_quadratic_n_elements: {}".format(fat_mesh_quadratic_n_elements))

variables.meshes["3DFatMesh_quadratic"] = {
  "nElements": fat_mesh_quadratic_n_elements,
  "nRanks": fat_mesh_n_ranks,
  "rankNos": variables.fat_global_rank_nos,
  "nodePositions": fat_mesh_node_positions_local,
  "inputMeshIsGlobal": False,
  "setHermiteDerivatives": False,
  "logKey": "3DFatMesh_quadratic"
}

debug = False
if False:
  print("settings 3DFatMesh: ")
  #print(str(variables.meshes["3Dmesh"]))
  with open("3Dmesh_{}".format(rank_no),"w") as f:
    f.write(str(variables.meshes["3Dmesh"]))

# coarsen quadratic meshes for elasticity
# -------------------------------------------
# At this point, we have the following meshes: 3Dmesh, 3Dmesh_quadratic, 3DFatMesh, 3DFatMesh_quadratic.
# (3Dmesh and 3Dmesh_quadratic) and also (3DFatMesh and 3DFatMesh_quadratic) share nodes, the linear elements are simply part of the quadratic elements and the partitioning is the same.
# Now, we coarsen the quadratic meshes which will be used for elasticity, because we do not need so much resolution.
# The coarser meshes will be 3Dmesh_elasticity_quadratic and 3DFatMesh_elasticity_quadratic.

# 3Dmesh_quadratic
# ----------------
n_elements_original_quadratic_mesh = variables.meshes["3Dmesh_quadratic"]["nElements"]    # local number of elements

import itertools

# determine number of elements
n_elements_aim_x = int(np.round(n_elements_original_quadratic_mesh[0] * (float)(variables.sampling_factor_elasticity_x)))
n_elements_aim_y = int(np.round(n_elements_original_quadratic_mesh[1] * (float)(variables.sampling_factor_elasticity_y)))
n_elements_aim_z = int(np.round(n_elements_original_quadratic_mesh[2] * (float)(variables.sampling_factor_elasticity_z)))

# crop to constraints, at least 1 element, maximum the available number
n_elements_aim_x = min(n_elements_original_quadratic_mesh[0], max(1, n_elements_aim_x))
n_elements_aim_y = min(n_elements_original_quadratic_mesh[1], max(1, n_elements_aim_y))
n_elements_aim_z = min(n_elements_original_quadratic_mesh[2], max(1, n_elements_aim_z))

# select nodes such that `n_elements_aim` quadratic 1D elements will be produced and they have their center point in the middle, e.g. not nodes [0,2,3,  5,6] (because element [0,2,3] is not good) but [0,2,4,  5,6] (element [0,2,4] is better)
node_indices_to_use_x = [int(np.round(x))*2 for x in np.linspace(0, n_elements_original_quadratic_mesh[0], n_elements_aim_x+1)]
node_indices_to_use_x = [(node_indices_to_use_x[i], int(0.5*(node_indices_to_use_x[i]+node_indices_to_use_x[i+1]))) for i in range(len(node_indices_to_use_x)-1)] + [(n_elements_original_quadratic_mesh[0]*2,)]
node_indices_to_use_x = list(itertools.chain(*node_indices_to_use_x))
#print("{}: selected node indices muscle mesh x: {}".format(rank_no,node_indices_to_use_x))

node_indices_to_use_y = [int(np.round(y))*2 for y in np.linspace(0, n_elements_original_quadratic_mesh[1], n_elements_aim_y+1)]
node_indices_to_use_y = [(node_indices_to_use_y[i], int(0.5*(node_indices_to_use_y[i]+node_indices_to_use_y[i+1]))) for i in range(len(node_indices_to_use_y)-1)] + [(n_elements_original_quadratic_mesh[1]*2,)]
node_indices_to_use_y = list(itertools.chain(*node_indices_to_use_y))

node_indices_to_use_z = [int(np.round(z))*2 for z in np.linspace(0, n_elements_original_quadratic_mesh[2], n_elements_aim_z+1)]
node_indices_to_use_z = [(node_indices_to_use_z[i], int(0.5*(node_indices_to_use_z[i]+node_indices_to_use_z[i+1]))) for i in range(len(node_indices_to_use_z)-1)] + [(n_elements_original_quadratic_mesh[2]*2,)]
node_indices_to_use_z = list(itertools.chain(*node_indices_to_use_z))

# determine number of nodes of the 3Dmesh_quadratic
n_points_local_original_quadratic_mesh_x = 2*n_elements_original_quadratic_mesh[0]
n_points_local_original_quadratic_mesh_y = 2*n_elements_original_quadratic_mesh[1]
n_points_local_original_quadratic_mesh_z = 2*n_elements_original_quadratic_mesh[2]

# if the own subdomain is at the (x+) boundary
if variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
  n_points_local_original_quadratic_mesh_x += 1

# if the own subdomain is at the (y+) boundary
if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
  n_points_local_original_quadratic_mesh_y += 1
  
# if the own subdomain is at the (z+) boundary
if variables.own_subdomain_coordinate_z == variables.n_subdomains_z - 1:
  n_points_local_original_quadratic_mesh_z += 1

# remove lasts indices if they are not part of the local nodes (ghost nodes)
if node_indices_to_use_x[-1] >= n_points_local_original_quadratic_mesh_x:
  #print(node_indices_to_use_x)
  node_indices_to_use_x = node_indices_to_use_x[:-1]
  #print("->",node_indices_to_use_x)
if node_indices_to_use_y[-1] >= n_points_local_original_quadratic_mesh_y:
  node_indices_to_use_y = node_indices_to_use_y[:-1]
if node_indices_to_use_z[-1] >= n_points_local_original_quadratic_mesh_z:
  node_indices_to_use_z = node_indices_to_use_z[:-1]

# get the node positions of the fine mesh
mesh_quadratic_node_positions_local = variables.meshes["3Dmesh_quadratic"]["nodePositions"]

if debug:
  print("{}: mesh {}".format(rank_no, "3Dmesh_quadratic"))
  print("{}: nElements old: {} x {} x {}".format(rank_no, n_elements_original_quadratic_mesh[0], n_elements_original_quadratic_mesh[1], n_elements_original_quadratic_mesh[2]))
  print("{}: nElements new: {} x {} x {}".format(rank_no, n_elements_aim_x, n_elements_aim_y, n_elements_aim_z))
  print("{}: n points old: {} x {} x {} = {}".format(rank_no, n_points_local_original_quadratic_mesh_x, n_points_local_original_quadratic_mesh_y, n_points_local_original_quadratic_mesh_z, n_points_local_original_quadratic_mesh_x*n_points_local_original_quadratic_mesh_y*n_points_local_original_quadratic_mesh_z))
  
  print("{}: n node positions old: {}".format(rank_no, len(mesh_quadratic_node_positions_local[1])))
  
node_positions_are_as_file_and_offset = isinstance(mesh_quadratic_node_positions_local[0], str)

# collect new node positions from fine node positions
# the indices of the nodes that should be used are in node_indices_to_use_{x,y,z}, given by the function `get_node_indices_to_use`
points_local_new_quadratic_mesh = []
for node_index_z in node_indices_to_use_z:
  for node_index_y in node_indices_to_use_y:
    for node_index_x in node_indices_to_use_x:
        
      index = node_index_z*n_points_local_original_quadratic_mesh_y*n_points_local_original_quadratic_mesh_x + node_index_y*n_points_local_original_quadratic_mesh_x + node_index_x
      
      #print("index {},{},{} = {}".format(node_index_x, node_index_y, node_index_z, index))
      
      # if fiber data was not loaded, the variable mesh_quadratic_node_positions_local is of type ['filename', [(filepos, length), ...]]
      if node_positions_are_as_file_and_offset:
        point = mesh_quadratic_node_positions_local[1][index]
      else:
        point = mesh_quadratic_node_positions_local[index]
        
      points_local_new_quadratic_mesh.append(point)
n_points_new = len(points_local_new_quadratic_mesh)
if debug:
  print("{}: n points new: {} x {} x {} = {}".format(rank_no,len(node_indices_to_use_x),len(node_indices_to_use_y),len(node_indices_to_use_z),n_points_new))

if node_positions_are_as_file_and_offset:
  points_local_new_quadratic_mesh = [mesh_quadratic_node_positions_local[0], points_local_new_quadratic_mesh]

# determine global number of nodes and elements for the "nPointsLocal", "nPointsGlobal" and "nElementsGlobal" options of the mesh
n_points_global_new_quadratic_mesh_x = 0
n_points_global_new_quadratic_mesh_y = 0
n_points_global_new_quadratic_mesh_z = 0
n_elements_global_new_quadratic_mesh_x = 0
n_elements_global_new_quadratic_mesh_y = 0
n_elements_global_new_quadratic_mesh_z = 0

# loop over subdomains in x direction
for subdomain_coordinate_x in range(variables.n_subdomains_x):
  n_elements_linear_x = n_sampled_points_in_subdomain_x(subdomain_coordinate_x)
  if subdomain_coordinate_x == variables.n_subdomains_x-1:
    n_elements_linear_x -= 1
  n_elements_quadratic_x = n_elements_linear_x // 2
  
  # determine n_elements_aim for the current subdomain in x direction
  n_quadratic_elements_aim_x = int(np.round(n_elements_quadratic_x * (float)(variables.sampling_factor_elasticity_x)))
  n_quadratic_elements_aim_x = min(n_elements_quadratic_x, max(1, n_quadratic_elements_aim_x))
  
  # determine number of nodes in the current subdomain
  n_nodes_x = 2*n_quadratic_elements_aim_x
  if subdomain_coordinate_x == variables.n_subdomains_x-1:
    n_nodes_x += 1
    
  # sum up global number of elements and nodes
  n_points_global_new_quadratic_mesh_x += n_nodes_x
  n_elements_global_new_quadratic_mesh_x += n_nodes_x//2
    
# loop over subdomains in y direction
for subdomain_coordinate_y in range(variables.n_subdomains_y):
  n_elements_linear_y = n_sampled_points_in_subdomain_y(subdomain_coordinate_y)
  if subdomain_coordinate_y == variables.n_subdomains_y-1:
    n_elements_linear_y -= 1
  n_elements_quadratic_y = n_elements_linear_y // 2
  
  # determine n_elements_aim for the current subdomain in y direction
  n_quadratic_elements_aim_y = int(np.round(n_elements_quadratic_y * (float)(variables.sampling_factor_elasticity_y)))
  n_quadratic_elements_aim_y = min(n_elements_quadratic_y, max(1, n_quadratic_elements_aim_y))
  
  # determine number of nodes in the current subdomain
  n_nodes_y = 2*n_quadratic_elements_aim_y
  if subdomain_coordinate_y == variables.n_subdomains_y-1:
    n_nodes_y += 1
    
  # sum up global number of elements and nodes
  n_points_global_new_quadratic_mesh_y += n_nodes_y
  n_elements_global_new_quadratic_mesh_y += n_nodes_y//2
    
# loop over subdomains in z direction
for subdomain_coordinate_z in range(variables.n_subdomains_z):
  n_elements_linear_z = n_sampled_points_in_subdomain_z(subdomain_coordinate_z)
  if subdomain_coordinate_z == variables.n_subdomains_z-1:
    n_elements_linear_z -= 1
  n_elements_quadratic_z = n_elements_linear_z // 2
  
  # determine n_elements_aim for the current subdomain in z direction
  n_quadratic_elements_aim_z = int(np.round(n_elements_quadratic_z * (float)(variables.sampling_factor_elasticity_z)))
  n_quadratic_elements_aim_z = min(n_elements_quadratic_z, max(1, n_quadratic_elements_aim_z))
  
  # determine number of nodes in the current subdomain
  n_nodes_z = 2*n_quadratic_elements_aim_z
  if subdomain_coordinate_z == variables.n_subdomains_z-1:
    n_nodes_z += 1
    
  # sum up global number of elements and nodes
  n_points_global_new_quadratic_mesh_z += n_nodes_z
  n_elements_global_new_quadratic_mesh_z += n_nodes_z//2
  
# set name of new mesh
elasticity_mesh_name = "3Dmesh_elasticity_quadratic"

variables.meshes[elasticity_mesh_name] = {
  "nElements":              [n_elements_aim_x, n_elements_aim_y, n_elements_aim_z],
  "nRanks":                 variables.meshes["3Dmesh_quadratic"]["nRanks"],
  "nodePositions":          points_local_new_quadratic_mesh,
  "inputMeshIsGlobal":      False,
  "setHermiteDerivatives":  False,
  "logKey":                 elasticity_mesh_name,
  
  # set information on how many nodes there are in the 3D mesh, this is not needed for the opendihu core but might be useful in some settings script or for debugging
  "nPointsLocal":           [len(node_indices_to_use_x), len(node_indices_to_use_y), len(node_indices_to_use_z)],
  "nPointsGlobal":          [n_points_global_new_quadratic_mesh_x, n_points_global_new_quadratic_mesh_y, n_points_global_new_quadratic_mesh_z],
  "nElementsGlobal":        [n_elements_global_new_quadratic_mesh_x, n_elements_global_new_quadratic_mesh_y, n_elements_global_new_quadratic_mesh_z]
}

if "rankNos" in variables.meshes["3Dmesh_quadratic"]:
  variables.meshes[elasticity_mesh_name]["rankNos"] = variables.meshes["3Dmesh_quadratic"]["rankNos"]

# 3DFatMesh_quadratic -> 3DFatMesh_elasticity_quadratic
# -------------------
# determine number of elements in fat elasticity mesh, x and z direction which is adjacent to the muscle mesh
# (in the analogous part for 3Dmesh, this is called n_elements_aim)
# on the top portion of the muscle mesh
if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
  n_elements_elasticity_fat_mesh_x = n_elements_aim_x
  
  # on the top right corner
  if variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
    n_elements_elasticity_fat_mesh_x = n_elements_aim_x + n_elements_aim_y

# on the right boundary but not the corner
elif variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
  n_elements_elasticity_fat_mesh_x = n_elements_aim_y
  
n_elements_elasticity_fat_mesh_z = n_elements_aim_z

# determine number of elements in y direction
n_elements_original_fat_mesh = variables.meshes["3DFatMesh_quadratic"]["nElements"]
n_elements_elasticity_fat_mesh_y = int(np.round(n_elements_original_fat_mesh[1] * (float)(variables.sampling_factor_elasticity_fat_y)))

# crop to constraints, at least 1 element, maximum the available number
n_elements_elasticity_fat_mesh_y = min(n_elements_original_fat_mesh[1], max(1, n_elements_elasticity_fat_mesh_y))

n_elements_elasticity_fat_mesh = [n_elements_elasticity_fat_mesh_x, n_elements_elasticity_fat_mesh_y, n_elements_elasticity_fat_mesh_z]

# in y direction, select which node indices from the original mesh should be used
# on the top portion of the muscle mesh
if variables.own_subdomain_coordinate_y == variables.n_subdomains_y - 1:
  node_indices_to_use_fat_x = list(node_indices_to_use_x)
  
  # on the top right corner
  if variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
    node_indices_to_use_fat_x += [n_points_local_original_quadratic_mesh_x-1+n_points_local_original_quadratic_mesh_y-1-i for i in reversed(node_indices_to_use_y)][1:]

# on the right boundary but not the corner
elif variables.own_subdomain_coordinate_x == variables.n_subdomains_x - 1:
  node_indices_to_use_fat_x = list(reversed(node_indices_to_use_y))
  
#print("{}: x-y-coord: ({},{})/({},{}), node indices in x direction, old x: {} y: {}, new: {}".format(rank_no, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.n_subdomains_x, variables.n_subdomains_y, \
#  node_indices_to_use_x,node_indices_to_use_y,node_indices_to_use_fat_x))


# in y direction, select nodes
node_indices_to_use_fat_y = [int(np.round(y))*2 for y in np.linspace(0, n_elements_original_fat_mesh[1], n_elements_elasticity_fat_mesh_y+1)]
node_indices_to_use_fat_y = [(node_indices_to_use_fat_y[i], int(0.5*(node_indices_to_use_fat_y[i]+node_indices_to_use_fat_y[i+1]))) for i in range(len(node_indices_to_use_fat_y)-1)] + [(n_elements_original_fat_mesh[1]*2,)]
node_indices_to_use_fat_y = list(itertools.chain(*node_indices_to_use_fat_y))

# determine all new positions of the fat mesh, using the node indices in x and z direction of the coarsened muscle mesh and node_indices_to_use_fat_y in y direction
points_local_new_quadratic_fat_mesh = []
for node_index_z in node_indices_to_use_z:
  for node_index_y in node_indices_to_use_fat_y:
    for node_index_x in node_indices_to_use_fat_x:
        
      index = node_index_z*fat_mesh_n_points[1]*fat_mesh_n_points[0] + node_index_y*fat_mesh_n_points[0] + node_index_x
      point = fat_mesh_node_positions_local[index]
        
      points_local_new_quadratic_fat_mesh.append(point)

points_local_new_fat_mesh = len(points_local_new_quadratic_fat_mesh)

variables.meshes["3DFatMesh_elasticity_quadratic"] = {
  "nElements":              n_elements_elasticity_fat_mesh,
  "nRanks":                 variables.meshes["3DFatMesh_quadratic"]["nRanks"],
  "nodePositions":          points_local_new_quadratic_fat_mesh,
  "inputMeshIsGlobal":      False,
  "setHermiteDerivatives":  False,
  "logKey":                 "3DFatMesh_elasticity_quadratic",
  
  # set information on how many nodes there are in the 3D mesh, this is not needed for the opendihu core but might be useful in some settings script or for debugging
  "nPointsLocal":           [len(node_indices_to_use_fat_x), len(node_indices_to_use_fat_y), len(node_indices_to_use_z)],
  "nPointsGlobal":          [n_points_global_new_quadratic_mesh_x+n_points_global_new_quadratic_mesh_y-1, 2*n_elements_elasticity_fat_mesh_y+1, n_points_global_new_quadratic_mesh_z],
  "nElementsGlobal":        [n_elements_global_new_quadratic_mesh_x+n_elements_global_new_quadratic_mesh_y, n_elements_elasticity_fat_mesh_y, n_elements_global_new_quadratic_mesh_z]
}

# output information about partitioning on rank 0
if rank_no == 0:      
  print("   sub-sampling 3D elasticity mesh with factors {}, {}, {} ".format(variables.sampling_factor_elasticity_x, variables.sampling_factor_elasticity_y, variables.sampling_factor_elasticity_z))
  
  n_points_global_elasticity_mesh = variables.meshes["3Dmesh_elasticity_quadratic"]["nPointsGlobal"]
  n_points_local_elasticity_mesh = variables.meshes["3Dmesh_elasticity_quadratic"]["nPointsLocal"]
  n_elements_global_elasticity_mesh = variables.meshes["3Dmesh_elasticity_quadratic"]["nElementsGlobal"]
  n_elements_local_elasticity_mesh = variables.meshes["3Dmesh_elasticity_quadratic"]["nElements"]
  n_points_global_elasticity_fat_mesh = variables.meshes["3DFatMesh_elasticity_quadratic"]["nPointsGlobal"]
  n_points_local_elasticity_fat_mesh = variables.meshes["3DFatMesh_elasticity_quadratic"]["nPointsLocal"]
  n_elements_global_elasticity_fat_mesh = variables.meshes["3DFatMesh_elasticity_quadratic"]["nElementsGlobal"]
  n_elements_local_elasticity_fat_mesh = variables.meshes["3DFatMesh_elasticity_quadratic"]["nElements"]
  print("   elasticity quadratic 3D meshes:")
  print("   muscle:             nodes global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
    n_points_global_elasticity_mesh[0], n_points_global_elasticity_mesh[1], n_points_global_elasticity_mesh[2], n_points_global_elasticity_mesh[0]*n_points_global_elasticity_mesh[1]*n_points_global_elasticity_mesh[2],
    n_points_local_elasticity_mesh[0], n_points_local_elasticity_mesh[1], n_points_local_elasticity_mesh[2], n_points_local_elasticity_mesh[0]*n_points_local_elasticity_mesh[1]*n_points_local_elasticity_mesh[2]))
  print("          quadratic elements global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
    n_elements_global_elasticity_mesh[0], n_elements_global_elasticity_mesh[1], n_elements_global_elasticity_mesh[2], n_elements_global_elasticity_mesh[0]*n_elements_global_elasticity_mesh[1]*n_elements_global_elasticity_mesh[2],
    n_elements_local_elasticity_fat_mesh[0], n_elements_local_elasticity_fat_mesh[1], n_elements_local_elasticity_fat_mesh[2], n_elements_local_elasticity_fat_mesh[0]*n_elements_local_elasticity_fat_mesh[1]*n_elements_local_elasticity_fat_mesh[2]))
  print("   fat and skin layer: nodes global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
    n_points_global_elasticity_fat_mesh[0], n_points_global_elasticity_fat_mesh[1], n_points_global_elasticity_fat_mesh[2], n_points_global_elasticity_fat_mesh[0]*n_points_global_elasticity_fat_mesh[1]*n_points_global_elasticity_fat_mesh[2],
    n_points_local_elasticity_fat_mesh[0], n_points_local_elasticity_fat_mesh[1], n_points_local_elasticity_fat_mesh[2], n_points_local_elasticity_fat_mesh[0]*n_points_local_elasticity_fat_mesh[1]*n_points_local_elasticity_fat_mesh[2]))
  print("          quadratic elements global: {} x {} x {} = {}, local: {} x {} x {} = {}".format( 
    n_elements_global_elasticity_fat_mesh[0], n_elements_global_elasticity_fat_mesh[1], n_elements_global_elasticity_fat_mesh[2], n_elements_global_elasticity_fat_mesh[0]*n_elements_global_elasticity_fat_mesh[1]*n_elements_global_elasticity_fat_mesh[2],
    n_elements_local_elasticity_fat_mesh[0], n_elements_local_elasticity_fat_mesh[1], n_elements_local_elasticity_fat_mesh[2], n_elements_local_elasticity_fat_mesh[0]*n_elements_local_elasticity_fat_mesh[1]*n_elements_local_elasticity_fat_mesh[2]))

#print("created meshes {}".format(list(variables.meshes.keys())))
# set output writers
# ---------------------------
# set output writer, this example doesn't use this, it defines the output writers directly in the settings file, but the code is left here in case we want to use it later when tidying up
if False:  
  variables.output_writer_fibers = []
  variables.output_writer_elasticity = []
  variables.output_writer_emg = []
  variables.output_writer_0D_states = []

  subfolder = ""
  if variables.paraview_output:
    if variables.adios_output:
      subfolder = "paraview/"
    variables.output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fixedFormat": False, "combineFiles": True})
    variables.output_writer_elasticity.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fixedFormat": False, "combineFiles": True})
    variables.output_writer_fibers.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True, "fixedFormat": False, "combineFiles": True})
    if variables.states_output:
      variables.output_writer_0D_states.append({"format": "Paraview", "outputInterval": 1, "filename": "out/" + subfolder + variables.scenario_name + "/0D_states", "binary": True, "fixedFormat": False, "combineFiles": True})

  if variables.adios_output:
    if variables.paraview_output:
      subfolder = "adios/"
    variables.output_writer_emg.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "useFrontBackBuffer": False, "combineNInstances": 1})
    variables.output_writer_elasticity.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "useFrontBackBuffer": False})
    variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "combineNInstances": variables.n_subdomains_xy, "useFrontBackBuffer": False})
    #variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + variables.scenario_name + "/fibers", "combineNInstances": 1, "useFrontBackBuffer": False}

  if variables.python_output:
    if variables.adios_output:
      subfolder = "python/"
    variables.output_writer_emg.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True})
    variables.output_writer_elasticity.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True})
    variables.output_writer_fibers.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True})

  if variables.exfile_output:
    if variables.adios_output:
      subfolder = "exfile/"
    variables.output_writer_emg.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg"})
    variables.output_writer_elasticity.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity"})
    variables.output_writer_fibers.append({"format": "Exfile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers"})

# set variable mappings for cellml model
if "hodgkin_huxley" in variables.cellml_file and "hodgkin_huxley-razumova" not in variables.cellml_file:
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
    ("connectorSlot","vm"):     "wal_environment/vS",                 # expose state 0 = Vm to the operator splitting
    ("connectorSlot", "stress"):"razumova/stress",                    # expose algebraic 12 = γ to the operator splitting
    ("connectorSlot", "lambda"):"razumova/L_S",                       # expose lambda
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
  
elif "hodgkin_huxley-razumova" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           "membrane/i_Stim",          # parameter 0 is I_stim
    ("parameter", 1):           "Razumova/l_hs",            # parameter 1 is fiber stretch λ
    ("connectorSlot", "vm"):    "membrane/V",               # expose Vm to the operator splitting
    ("connectorSlot", "stress"):"Razumova/activestress",
    ("connectorSlot", "alpha"): "Razumova/activation",      # expose activation .
  }
  variables.parameters_initial_values = [0, 1]
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

else:
  print("\033[0;31mCellML file {} has no mappings implemented in helper.py\033[0m".format(variables.cellml_file))
  quit()

# load MU distribution and firing times
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# ---------------------
# callback functions
def get_motor_unit_no(compartment_no):
  return compartment_no

def compartment_gets_stimulated(compartment_no, frequency, current_time):
  """
  determine if compartment compartment_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(compartment_no)*alpha)
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(variables.firing_times,0)
  
  #if variables.firing_times[index % n_firing_times, mu_no] == 1:
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), variables.firing_times[index % n_firing_times, mu_no], "true" if variables.firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return variables.firing_times[index % n_firing_times, mu_no] == 1
  
# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, compartment_no):
  
  # determine if fiber gets stimulated at the current time
  is_compartment_gets_stimulated = compartment_gets_stimulated(compartment_no, variables.stimulation_frequency, current_time)

  if is_compartment_gets_stimulated:  
    
    n_nodes_x = variables.n_points_3D_mesh_global_x
    n_nodes_y = variables.n_points_3D_mesh_global_y
    n_nodes_z = variables.n_points_3D_mesh_global_z
    z_index_center = (int)(n_nodes_z/2)
    y_index_center = (int)(n_nodes_y/2)
    x_index_center = (int)(n_nodes_x/2)
    
    for k in range(n_nodes_z):
      if z_index_center-1 <= k <= z_index_center+1:  # use only 3 nodes in z direction from center
        for j in range(n_nodes_y):
          #if y_index_center-1 <= j <= y_index_center+1:  # use only 3 nodes in y direction from center
          if True:                                        # use all nodes
            for i in range(n_nodes_x):
              #if x_index_center-1 <= i <= x_index_center+1:  # use only 3 nodes in x direction from center
              if True:                                        # use all nodes
                
                key = ((i,j,k),0,0)        # key: ((x,y,z),nodal_dof_index,state_no)
                states[key] = variables.vm_value_stimulated
                #print("set states at ({},{},{}) to 40".format(i,j,k))

    #print("states: {}".format(states))
    #print("n_nodes: ({},{},{})".format(n_nodes_x, n_nodes_y, n_nodes_z))
    #print("n_nodes_global: {}, time_step_no: {}, current_time: {}, compartment_no: {}".format(n_nodes_global, time_step_no, current_time, compartment_no))
    #wait = input("Press any key to continue...")

def set_dummy_stress(n_dofs_global, time_step_no, current_time, values, global_natural_dofs, compartment_no):
  pass
  
# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
  print("\nDebugging output about compartment firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(variables.firing_times,0)
  for compartment_no_index in range(variables.n_compartments):
    if compartment_no_index % 100 == 0:
      t_algebraic = timeit.default_timer()
      if t_algebraic - t_start > 100:
        print("Note: break after {}/{} compartments ({:.0f}%) because it already took {:.3f}s".format(compartment_no_index,variables.n_compartments,100.0*compartment_no_index/(variables.n_compartments-1.),t_algebraic - t_start))
        break
    
    first_stimulation = None
    for current_time in np.linspace(0,1./variables.stimulation_frequency*n_firing_times,n_firing_times):
      if compartment_gets_stimulated(compartment_no_index, variables.stimulation_frequency, current_time):
        first_stimulation = current_time
        break
    mu_no = get_motor_unit_no(compartment_no_index)
    first_stimulation_info.append([compartment_no_index,mu_no,first_stimulation])
  
  first_stimulation_info.sort(key=lambda x: 1e6+1e-6*x[1]+1e-12*x[0] if x[2] is None else x[2]+1e-6*x[1]+1e-12*x[0])
  
  print("First stimulation times")
  print("    Time  MU compartments")
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
  
####################################
# set Dirichlet BC for the flow problem

n_points_3D_mesh_linear_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_linear_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_linear_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_linear_global = n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y*n_points_3D_mesh_linear_global_z

n_points_3D_mesh_quadratic_global_x = 2*n_points_3D_mesh_linear_global_x - 1
n_points_3D_mesh_quadratic_global_y = 2*n_points_3D_mesh_linear_global_y - 1
n_points_3D_mesh_quadratic_global_z = 2*n_points_3D_mesh_linear_global_z - 1
 
# set Dirichlet BC values for bottom nodes to 0 and for top nodes to 1
variables.potential_flow_dirichlet_bc = {}
for i in range(n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y):
  variables.potential_flow_dirichlet_bc[i] = 0.0
  variables.potential_flow_dirichlet_bc[(n_points_3D_mesh_linear_global_z-1)*n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y + i] = 1.0

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y
  
# set boundary conditions for the elasticity
# Note, we have a composite mesh, consisting of 3Dmesh_elasticity_quadratic and 3DFatMesh_elasticity_quadratic and this composite mesh has a numbering that goes over all dofs.
# The following works because we index the first sub mesh and the first mesh of a composite mesh always has all own dofs with their normal no.s. (The 2nd mesh has the shared dofs to the first mesh removed in the numbering, i.e. they are not counted twice).
[mx, my, mz] = variables.meshes["3Dmesh_elasticity_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_elasticity_quadratic"]["nElementsGlobal"]

# print all mesh information
if False:
  for mesh_name in ["3Dmesh_elasticity_quadratic","3DFatMesh_elasticity_quadratic"]:
    print("")
    print("{}: mesh {}".format(rank_no, mesh_name))
    print("    nPointsLocal: {}".format(variables.meshes[mesh_name]["nPointsLocal"]))
    print("    nPointsGlobal: {}".format(variables.meshes[mesh_name]["nPointsGlobal"]))
    print("    nElementsGlobal: {}".format(variables.meshes[mesh_name]["nElementsGlobal"]))
    print("    inputMeshIsGlobal: {}".format(variables.meshes[mesh_name]["inputMeshIsGlobal"]))
    print("    nElements: {}".format(variables.meshes[mesh_name]["nElements"]))
    print("    nRanks: {}".format(variables.meshes[mesh_name]["nRanks"]))
    print("    len(nodePositions): {}".format(len(variables.meshes[mesh_name]["nodePositions"][1])))

# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top
variables.elasticity_dirichlet_bc = {}
for j in range(my):
  for i in range(mx):
    variables.elasticity_dirichlet_bc[(mz-1)*mx*my + j*mx + i] = [None,None,0.0,None,None,None]
  
# fix edge
#for i in range(mx):
#  variables.elasticity_dirichlet_bc[(mz-1)*mx*my + 0*mx + i] = [0.0,None,0.0,None,None,None]
  
# fix corner completely
variables.elasticity_dirichlet_bc[(mz-1)*mx*my + 0] = [0.0,0.0,0.0,None,None,None]

# Neumann BC at bottom nodes, traction downwards
variables.elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": variables.bottom_traction, "face": "2-"} for j in range(ny) for i in range(nx)]
#variables.elasticity_neumann_bc = []

####################################
# compute relative factors fr for compartments
def compute_compartment_relative_factors(mesh_node_positions, n_mesh_points_xy, n_mesh_points_z, fiber_data, motor_units):
  """
  Compute the relative factors, f_r, that are needed in the multidomain formulation as a weighting for compartments.
  Result is relative_factors[motor_unit_no][node_no] for the 3D mesh.
  :param mesh_node_positions:  list of (x,y,z) values, global node positions of the 3D mesh
  :param fiber_data: list of fibers, each fiber is a list of points, i.e. point = fiber_data[xy_index][z_index]
  :param motor_units: a list of dicts, settings for the motor units, [{"fiber_no": 0, "standard_deviation": 0.5, "maximum": 1}]
  """
  
  # list of fibers, fiber = list of points, point = list with 3 coordinate entries
  n_compartments = len(motor_units)
  n_points_fiber = len(fiber_data[0])

  # create relative factors for compartments
  #if rank_no == 0:
  #  print("determine relative factors for {} motor units:\n{}".format(n_compartments, motor_units))

  # determine approximate diameter of muscle at every point is z direction
  diameters = []
    
  # loop over points in z direction
  for z_index_mesh in range(n_mesh_points_z):
    
    z_index_fiber = int(z_index_mesh / (float)(n_mesh_points_z) * n_points_fiber)
    
    # get point on first and last fiber
    point0 = np.array(fiber_data[0][z_index_fiber])
    point4 = np.array(fiber_data[(variables.n_fibers_x-1)//2][z_index_fiber])
    point1 = np.array(fiber_data[variables.n_fibers_x-1][z_index_fiber])
    point2 = np.array(fiber_data[-variables.n_fibers_x][z_index_fiber])
    point5 = np.array(fiber_data[(-variables.n_fibers_x)//2][z_index_fiber])
    point3 = np.array(fiber_data[-1][z_index_fiber])
    
    # their distance is an approximation for the diameter
    distance01 = np.linalg.norm(point0 - point1)
    distance02 = np.linalg.norm(point0 - point2)
    distance03 = np.linalg.norm(point0 - point3)
    distance04 = np.linalg.norm(point0 - point4)
    distance05 = np.linalg.norm(point0 - point5)
    distance12 = np.linalg.norm(point1 - point2)
    distance13 = np.linalg.norm(point1 - point3)
    distance14 = np.linalg.norm(point1 - point4)
    distance15 = np.linalg.norm(point1 - point5)
    distance23 = np.linalg.norm(point2 - point3)
    distance24 = np.linalg.norm(point2 - point4)
    distance25 = np.linalg.norm(point2 - point5)
    distance34 = np.linalg.norm(point3 - point4)
    distance35 = np.linalg.norm(point3 - point5)
    distance45 = np.linalg.norm(point4 - point5)
    distance = max(distance01,distance02,distance03,distance04,distance05,distance12,distance13,distance14,distance15,distance23,distance24,distance25,distance34,distance35,distance45)
    diameters.append(distance)

  #print("diameters: {}".format(diameters))

  # create data structure with 0
  relative_factors = np.zeros((n_compartments, len(mesh_node_positions)))   # each row is one compartment

  # loop over nodes of mesh
  for node_no,node_position in enumerate(mesh_node_positions):
    node_position = np.array(node_position)
    
    z_index_mesh = int((float)(node_no) / n_mesh_points_xy)
    z_index_fiber = int(z_index_mesh / (float)(n_mesh_points_z) * n_points_fiber)
    
    # loop over motor units
    for motor_unit_no,motor_unit in enumerate(motor_units):
      
      # find point on fiber that is closest to current node
      fiber_no = motor_unit["fiber_no"]
      if fiber_no >= len(fiber_data):
        new_fiber_no = fiber_no % len(fiber_data)
        if node_no == 0:
          print("\033[0;31mError with motor unit {} around fiber {}, only {} fibers available, now using fiber {} % {} = {} instead.\033[0m".format(motor_unit_no, fiber_no, len(fiber_data), fiber_no, len(fiber_data), new_fiber_no))
        fiber_no = new_fiber_no
      
      min_distance = None
      search_range = int(1 / (float)(n_mesh_points_z) * n_points_fiber)
      search_range = max(10,search_range)
      z_start = max(0,z_index_fiber - search_range)
      z_end = min(n_points_fiber, z_index_fiber + search_range)
      
      #print("node_position: {}, z_index_fiber: {}, fiber at z index: {}, fiber: {}".format(node_position, z_index_fiber, fiber_data[fiber_no][z_index_fiber], fiber_data[fiber_no][z_start:z_end]))
      #print("search_range: {}".format(search_range))
      
      for k,fiber_point in enumerate(fiber_data[fiber_no][z_start:z_end]):
        d = np.array(fiber_point) - node_position
        distance = np.inner(d,d)
        if min_distance is None or distance < min_distance:
          min_distance = distance
          #print("node_position {}, fiber_point {}, d={}, |d|={}".format(node_position, fiber_point, d, np.sqrt(distance)))
      
      distance = np.sqrt(min_distance)
      
      # compute value as gaussian with given standard_deviation and maximum
      standard_deviation = motor_unit["standard_deviation"]*diameters[z_index_mesh]
      gaussian = scipy.stats.norm(loc = 0., scale = standard_deviation)
      value = gaussian.pdf(distance)*standard_deviation*np.sqrt(2*np.pi)*motor_unit["maximum"]
      relative_factors[motor_unit_no][node_no] += value
      
      #print("motor unit {}, fiber {}, distance {}, value {}".format(motor_unit_no, fiber_no, distance, value))

  return relative_factors

####################################
# load relative factors for motor units

# determine relative factor fields fr(x) for compartments
if not os.path.exists(variables.relative_factors_file):

  # the file does not yet exist, create it on rank 0
  if rank_no == 0: 
    
    mesh_node_positions = variables.meshes["3Dmesh"]["globalNodePositions"]
    n_points_global = variables.meshes["3Dmesh"]["nPointsGlobal"]
    n_mesh_points_xy = n_points_global[0]*n_points_global[1]
    n_mesh_points_z = n_points_global[2]
    
    print("Computing the relative MU factors, f_r, for {} motor units and {} mesh nodes, {} fibers. This may take a while ...".format(len(variables.motor_units), len(mesh_node_positions), len(variables.fibers)))
    variables.relative_factors = compute_compartment_relative_factors(mesh_node_positions, n_mesh_points_xy, n_mesh_points_z, variables.fibers, variables.motor_units)
    if rank_no == 0:
      print("Save relative factors to file \"{}\".".format(variables.relative_factors_file))
      with open(variables.relative_factors_file, "wb") as f:
        pickle.dump(variables.relative_factors, f)
  else:
    # wait until file is created on rank 0
    while not os.path.exists(variables.relative_factors_file):
      time.sleep(1)

if os.path.exists(variables.relative_factors_file):
  with open(variables.relative_factors_file, "rb") as f:
    if rank_no == 0:
      print("Load relative factors, f_r, from file \"{}\"".format(variables.relative_factors_file))
    variables.relative_factors = pickle.load(f, encoding='latin1')
else:
  print("\033[0;31mError: Could not load relative factors file \"{}\"\033[0m".format(variables.relative_factors_file))
  quit()

# debugging output
if rank_no == 0 and not variables.disable_firing_output:
  for i,factors_list in enumerate(variables.relative_factors.tolist()):
    print("MU {}, maximum fr: {}".format(i,max(factors_list)))

    
