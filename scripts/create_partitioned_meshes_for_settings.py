# Helper script to create the mesh partitioning that is needed for config Dict for opendihu
#
# The main function in this script is create_partitioned_meshes_for_settings, which creates the actual partitioning. 
# It returns some information about the partitioning, e.g. the local node positions and elements. 
# Furthermore it sets some more variables in the `variables` module, e.g. `variables.n_fibers_per_subdomain_x`.
# This script also contains functions, that give more information about the local portion of the partitioning, i.e. the local subdomain.
# These functions need the variables that are set by create_partitioned_meshes_for_settings.

import numpy as np
import pickle
import sys
import struct
sys.path.insert(0, '..')
import variables    # file variables.py

# define helper functions for fiber numbering

# number of fibers that are handled inside the subdomain x
def n_fibers_in_subdomain_x(subdomain_coordinate_x):
  a1 = variables.n_fibers_x - variables.n_subdomains_x*variables.n_fibers_per_subdomain_x              # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_x - a1                                               # number of subdomains with low number of fibers
  if subdomain_coordinate_x < a1:
    return variables.n_fibers_per_subdomain_x+1      # high number of fibers
  else:
    return variables.n_fibers_per_subdomain_x    # low number of fibers
  
# number of fibers that are handled inside the subdomain y
def n_fibers_in_subdomain_y(subdomain_coordinate_y):
  a1 = variables.n_fibers_y - variables.n_subdomains_y*variables.n_fibers_per_subdomain_y              # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_y - a1                                               # number of subdomains with low number of fibers
  if subdomain_coordinate_y < a1:
    return variables.n_fibers_per_subdomain_y+1     # high number of fibers
  else:
    return variables.n_fibers_per_subdomain_y     # low number of fibers

def n_fibers_in_previous_subdomains_y(subdomain_coordinate_y):
  # number of fibers handled in previous subdomains in y direction
  a1 = variables.n_fibers_y - variables.n_subdomains_y*variables.n_fibers_per_subdomain_y              # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_y - a1                                               # number of subdomains with low number of fibers
  
  if subdomain_coordinate_y < a1:
    return subdomain_coordinate_y * (variables.n_fibers_per_subdomain_y+1)
  else:
    return a1 * (variables.n_fibers_per_subdomain_y+1) + (subdomain_coordinate_y-a1) * variables.n_fibers_per_subdomain_y
    
def n_fibers_in_previous_subdomains_x(subdomain_coordinate_x):
  # number of fibers handled in previous subdomains in x direction
  a1 = variables.n_fibers_x - variables.n_subdomains_x*variables.n_fibers_per_subdomain_x              # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_x - a1                                               # number of subdomains with low number of fibers
  
  if subdomain_coordinate_x < a1:
    return subdomain_coordinate_x * (variables.n_fibers_per_subdomain_x+1)
  else:
    return a1 * (variables.n_fibers_per_subdomain_x+1) + (subdomain_coordinate_x-a1) * variables.n_fibers_per_subdomain_x

# global fiber no, from subdomain coordinate and coordinate inside the subdomain
def fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y):
  return (n_fibers_in_previous_subdomains_y(subdomain_coordinate_y) + fiber_in_subdomain_coordinate_y)*variables.n_fibers_x \
    + n_fibers_in_previous_subdomains_x(subdomain_coordinate_x) + fiber_in_subdomain_coordinate_x

# number of points that are handled inside the subdomain z (without ghost points)
def n_points_in_subdomain_z(subdomain_coordinate_z):
  a1 = variables.n_points_whole_fiber - variables.n_subdomains_z*variables.n_points_per_subdomain_z              # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_z - a1                                               # number of subdomains with low number of fibers
  if subdomain_coordinate_z < a1:
    return variables.n_points_per_subdomain_z+1     # high number of points
  else:
    return variables.n_points_per_subdomain_z     # low number of points
  
def n_points_in_previous_subdomains_z(subdomain_coordinate_z):
  # number of points handled in previous subdomains in z direction
  a1 = variables.n_points_whole_fiber - variables.n_subdomains_z*variables.n_points_per_subdomain_z              # number of subdomains with high number of points
  a2 = variables.n_subdomains_z - a1                                               # number of subdomains with low number of points
  
  if subdomain_coordinate_z < a1:
    return subdomain_coordinate_z * (variables.n_points_per_subdomain_z+1)
  else:
    return a1 * (variables.n_points_per_subdomain_z+1) + (subdomain_coordinate_z-a1) * variables.n_points_per_subdomain_z

# number of points in 3D mesh in the subdomain, in x direction
def n_sampled_points_in_subdomain_x(subdomain_coordinate_x):
  result = (int)(np.ceil(n_fibers_in_subdomain_x(subdomain_coordinate_x) / variables.sampling_stride_x))
  if subdomain_coordinate_x == variables.n_subdomains_x-1 and (n_fibers_in_subdomain_x(subdomain_coordinate_x)-1) % variables.sampling_stride_x != 0:
    result += 1
  return result

# number of points in 3D mesh in the subdomain, in y direction
def n_sampled_points_in_subdomain_y(subdomain_coordinate_y):
  result = (int)(np.ceil(n_fibers_in_subdomain_y(subdomain_coordinate_y) / variables.sampling_stride_y))
  if subdomain_coordinate_y == variables.n_subdomains_y-1 and (n_fibers_in_subdomain_y(subdomain_coordinate_y)-1) % variables.sampling_stride_y != 0:
    result += 1
  return result

# number of points in 3D mesh in the subdomain, in z direction
def n_sampled_points_in_subdomain_z(subdomain_coordinate_z):
  result = (int)(np.ceil(n_points_in_subdomain_z(subdomain_coordinate_z) / variables.sampling_stride_z))
  if subdomain_coordinate_z == variables.n_subdomains_z-1 and (n_points_in_subdomain_z(subdomain_coordinate_z)-1) % variables.sampling_stride_z != 0:
    result += 1
  return result

# this is the main function that creates the meshes and partitioning
def create_partitioned_meshes_for_settings(n_subdomains_x, n_subdomains_y, n_subdomains_z, fiber_file, load_fiber_data, sampling_stride_x, sampling_stride_y, sampling_stride_z, quadratic_3d_mesh=False):
  """
  Parse the binary fiber geometry file and creates 1D fiber meshes and a 3D mesh (linear or quadratic elements for the 3D mesh)
  It create the `meshes` Dict that can directly be used in the Python config like: 
  config = {
    "Meshes": meshes,
    ...
  }
  
  One 3D mesh with key either "3Dmesh" or "3Dmesh_quadratic" will be created, and multiple 1D fiber meshes with keys "MeshFiber_0", "MeshFiber_1", etc.
  The 3D mesh is aligned with the 1D meshes.
  
  The mesh specifications have the following form:
  A 1D mesh:
  meshes["MeshFiber_0"] = {
    "nElements":             n_fiber_elements_on_subdomain,    # the number of 1D elements of this fiber on the own processes' subdomain
    "nodePositions":         fiber_node_positions,             # all node positions of this fiber on the own processes' subdomain
    "inputMeshIsGlobal":     False,                            # `False` means that only locale nodes and elements are provided
    "nRanks":                [variables.n_subdomains_z],       # the number of ranks in z-direction over which this fiber mesh is spread
    "setHermiteDerivatives": False
  }
  
  The 3D mesh:
  meshes["3Dmesh"] = {
    "nElements":             variables.n_elements_3D_mesh,     # the number of elements on the own process
    "nRanks":                [n_subdomains_x, n_subdomains_y, n_subdomains_z],   # the number of subdomains = processes of the 3D mesh
    "nodePositions":         node_positions_3d_mesh,           # all local node positions
    "inputMeshIsGlobal":     False,                            # `False` means that only locale nodes and elements are provided
    "setHermiteDerivatives": False,
    "logKey": "3Dmesh"
  }
  
  The quadratic 3D mesh: same as 3D mesh, but under meshes["3Dmesh_quadratic"], with same number of elements but with quadratic elements with 27 nodes instead of 8, 
  the missing nodes are interpolated.
  
  This means that all quantities like, e.g. the local 3D node positions can be retrieved by meshes["3Dmesh"]["nodePositions"] and so on.
  
  :param n_subdomains_x: number of subdomains in x direction.
  :param n_subdomains_y: number of subdomains in y direction
  :param n_subdomains_z: number of subdomains in z direction, i.e. the number of subdivisions per fiber. 
                         The z direction is the direction of the fibers along the muscle.
                         the product n_subdomains_x*n_subdomains_y*n_subdomains_z has to equal the number of processes
  :param fiber_file: file name of the fiber file "*.bin" which contains fiber geometry in 3D space. The file has a pickle format.
  :param load_fiber_data: If the actual geometry data should be read from the file and the node positions placed in the generated meshes Dict. 
                          If this is true, the node positions will be in the mesh config.
                          If this is false, it inserts the filename and the position/offset in the file where the geometry data is present.
                          This information is later used by the C++ code to read the file in parallel. This is necessary for very large runs,
                          where the file has to be parsed in parallel and not already here in the python script.
  :param sampling_stride_x: Grid point stride in x direction to make the 3D mesh coarser than the grid points of the 1D fibers.
                            E.g. 2 means there will be 2 fibers per 3D element in x-direction on average.
                            (Or 3 if you want, two of them are on the edges of the 3D elements, one is along the center)
  :param sampling_stride_y: Grid point stride in y direction to make the 3D mesh coarser than the grid points of the 1D fibers.
  :param sampling_stride_z: Grid point stride in z direction to make the 3D mesh coarser than the grid points of the 1D fibers.
                            E.g. 20 means that there will be 20 mesh nodes of the 1D mesh per 3D element. 
                            The z direction is the direction of the fibers along the muscle.
  :param quadratic_3d_mesh: Whether to create the linear or the quadratic mesh (under "3Dmesh" or "3Dmesh_quadratic"), only one of both
                            will be created.
  
  :return: a list of the following entries:
  [meshes, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, n_fibers_x, n_fibers_y, n_points_whole_fiber],
  
  - meshes contains the Meshes as python dicts, as explained above
  - where own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z are the cartesian coordinates of
    the own subdomain in the range {0,...,n_subdomains_x-1}, {0,...,n_subdomains_y-1}, {0,...,n_subdomains_z-1}.
  - n_fibers_x, n_fibers_y are the number of fibers in x and y direction that were parsed from the fiber file
  - n_points_whole_fiber: the number of grid points in every fiber (the number is the same for every fiber)
  
  
  """

  # store parameters under `variables` module
  variables.n_subdomains_x = n_subdomains_x
  variables.n_subdomains_y = n_subdomains_y
  variables.n_subdomains_z = n_subdomains_z
  variables.sampling_stride_x = sampling_stride_x
  variables.sampling_stride_y = sampling_stride_y
  variables.sampling_stride_z = sampling_stride_z

  # parse arguments concerning number of MPI ranks and own MPI rank
  # these values are passed from opendihu
  rank_no = (int)(sys.argv[-2])
  n_ranks = (int)(sys.argv[-1])

  variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
  own_subdomain_coordinate_x = rank_no % variables.n_subdomains_x
  own_subdomain_coordinate_y = (int)(rank_no / variables.n_subdomains_x) % variables.n_subdomains_y
  own_subdomain_coordinate_z = (int)(rank_no / variables.n_subdomains_xy)
  #print("rank: {}/{}".format(rank_no,n_ranks))

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

  # parse parameters in the file
  parameters = []
  for i in range(int(header_length/4.) - 1):
    double_raw = fiber_file_handle.read(4)
    value = struct.unpack('i', double_raw)[0]
    parameters.append(value)
    
  variables.n_fibers_total = parameters[0]
  variables.n_fibers_x = (int)(np.round(np.sqrt(variables.n_fibers_total)))
  variables.n_fibers_y = variables.n_fibers_x
  variables.n_points_whole_fiber = parameters[1]

  if rank_no == 0:
    print("n fibers:              {} ({} x {})".format(variables.n_fibers_total, variables.n_fibers_x, variables.n_fibers_y))
    print("n points per fiber:    {}".format(variables.n_points_whole_fiber))
      
  # parse whole fiber file, only if enabled
  if load_fiber_data:
    variables.fibers = []
    for fiber_index in range(variables.n_fibers_total):
      fiber = []
      for point_no in range(variables.n_points_whole_fiber):
        point = []
        for i in range(3):
          double_raw = fiber_file_handle.read(8)
          value = struct.unpack('d', double_raw)[0]
          point.append(value)
        fiber.append(point)
      variables.fibers.append(fiber)
           
  # compute partitioning
  if rank_no == 0:
    if n_ranks != variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z:
      print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
      quit()
    
  variables.n_fibers_per_subdomain_x = (int)(variables.n_fibers_x / variables.n_subdomains_x)
  variables.n_fibers_per_subdomain_y = (int)(variables.n_fibers_y / variables.n_subdomains_y)
  variables.n_points_per_subdomain_z = (int)(variables.n_points_whole_fiber / variables.n_subdomains_z)

  #####################
  # define 3D mesh

  n_sampled_points_in_own_subdomain_x = n_sampled_points_in_subdomain_x(own_subdomain_coordinate_x)
  n_sampled_points_in_own_subdomain_y = n_sampled_points_in_subdomain_y(own_subdomain_coordinate_y)
  n_sampled_points_in_own_subdomain_z = n_sampled_points_in_subdomain_z(own_subdomain_coordinate_z)
  
  # determine node positions of the 3D mesh
  node_positions_3d_mesh = []
  if not load_fiber_data:
    node_positions_3d_mesh = [fiber_file, []]  # first component: fiber file name, second component: data chunks, chunk = (offset, number of points)

  # range of points in z direction
  variables.z_point_index_start = n_points_in_previous_subdomains_z(own_subdomain_coordinate_z)
  variables.z_point_index_end = variables.z_point_index_start + n_points_in_subdomain_z(own_subdomain_coordinate_z)

  # loop over z point indices
  for k in range(n_sampled_points_in_own_subdomain_z):
    z_point_index = variables.z_point_index_start + k*variables.sampling_stride_z
    
    if own_subdomain_coordinate_z == variables.n_subdomains_z-1 and k == n_sampled_points_in_own_subdomain_z-1:
      z_point_index = variables.z_point_index_end-1
      
    #print("{}: sampling_stride_z: {}, k: {}, z: {}/{}".format(rank_no, variables.sampling_stride_z, k, z_point_index, variables.z_point_index_end))
    
    # loop over variables.fibers for own rank
    # loop over fiber in y-direction
    for j in range(n_sampled_points_in_own_subdomain_y):
      fiber_in_subdomain_coordinate_y = j*variables.sampling_stride_y
      
      # on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
      if own_subdomain_coordinate_y == variables.n_subdomains_y-1 and j == n_sampled_points_in_own_subdomain_y-1:
        fiber_in_subdomain_coordinate_y = n_fibers_in_subdomain_y(own_subdomain_coordinate_y)-1
      
      #print("{}: sampling_stride_y: {}, j: {}, y: {}/{}".format(rank_no, variables.sampling_stride_y, j, fiber_in_subdomain_coordinate_y, n_fibers_in_subdomain_y(own_subdomain_coordinate_y)))
        
      # loop over fiber in x-direction
      for i in range(n_sampled_points_in_own_subdomain_x):
        fiber_in_subdomain_coordinate_x = i*variables.sampling_stride_x
        
        # on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
        if own_subdomain_coordinate_x == variables.n_subdomains_x-1 and i == n_sampled_points_in_own_subdomain_x-1:
          fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(own_subdomain_coordinate_x)-1
        
        #print("{}: sampling_stride_x: {}, i: {}, x: {}/{}".format(rank_no, variables.sampling_stride_x, i, fiber_in_subdomain_coordinate_x, n_fibers_in_subdomain_x(own_subdomain_coordinate_x)))
        
        # get fiber no
        fiber_index = fiber_no(own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
        
        # read point from fiber file
        memory_size_fiber = variables.n_points_whole_fiber*3*8
        offset = 32 + header_length + fiber_index*memory_size_fiber + z_point_index*3*8
        
        if load_fiber_data:
          
          fiber_file_handle.seek(offset)
        
          point = []
          for j in range(3):
            double_raw = fiber_file_handle.read(8)
            value = struct.unpack('d', double_raw)[0]
            point.append(value)
        
          reference_point = variables.fibers[fiber_index][z_point_index]
          
          difference = np.linalg.norm(np.array(reference_point) - np.array(point))
          if difference > 1e-3:
            print("Error, point does not match: reference_point: ", reference_point, ", point: ", point)
            quit()
          node_positions_3d_mesh.append(point)
          
        else:
          node_positions_3d_mesh[1].append((offset, 1))   # command to read 1 point from offset
        #print("{}: {}".format(rank_no, point))
       
       
  # set local number of elements for the 3D mesh
  variables.n_elements_3D_mesh = [
      n_sampled_points_in_own_subdomain_x,
      n_sampled_points_in_own_subdomain_y, 
      n_sampled_points_in_own_subdomain_z
    ]
  
  # border subdomains have one element less than fibers
  if own_subdomain_coordinate_x == variables.n_subdomains_x-1:
    variables.n_elements_3D_mesh[0] -= 1
  if own_subdomain_coordinate_y == variables.n_subdomains_y-1:
    variables.n_elements_3D_mesh[1] -= 1
  if own_subdomain_coordinate_z == variables.n_subdomains_z-1:
    variables.n_elements_3D_mesh[2] -= 1

  print("len: {}".format(len(node_positions_3d_mesh)))
  print("n_elements_3D_mesh: {}".format(variables.n_elements_3D_mesh))
  
  # set the entry for the config
  meshes = {}
  if not quadratic_3d_mesh:
    meshes["3Dmesh"] = {
      "nElements": variables.n_elements_3D_mesh,
      "nRanks": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z],
      "nodePositions": node_positions_3d_mesh,
      "inputMeshIsGlobal": False,
      "setHermiteDerivatives": False,
      "logKey": "3Dmesh"
    }
  
  n_points_3D_mesh_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
  n_points_3D_mesh_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
  n_points_3D_mesh_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
  
  # create quadratic 3D mesh
  if quadratic_3d_mesh:
    node_positions_3d_mesh_quadratic = []
    
    # loop over existing linear 3D elements
    for k in range(n_sampled_points_in_own_subdomain_z):
      for j in range(n_sampled_points_in_own_subdomain_y):
        #   +--+
        #  /  /|
        # +--+ *1
        # |  |*
        # +--*0
        for i in range(n_sampled_points_in_own_subdomain_x):
          point0 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i]
          node_positions_3d_mesh_quadratic.append(point0)
          
          if i < n_sampled_points_in_own_subdomain_x-1:
            point1 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i+1]
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point1))/2.0))
            
        #   +--+
        #  /| /|
        # +-3* 1 (bottom center)
        # | *|/
        # 2**0
        for i in range(n_sampled_points_in_own_subdomain_x):
          point0 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i]
          point1 = None
          point2 = None
          point3 = None
          
          if i < n_sampled_points_in_own_subdomain_x-1:
            point1 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i+1]
          
          if j < n_sampled_points_in_own_subdomain_y-1:
            point2 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + (j+1)*n_sampled_points_in_own_subdomain_x + i]
            
          if i < n_sampled_points_in_own_subdomain_x-1 and j < n_sampled_points_in_own_subdomain_y-1:
            point3 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + (j+1)*n_sampled_points_in_own_subdomain_x + i+1]
          
          if point0 is not None and point2 is not None:
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point2))/2.0))
            
          if point0 is not None and point1 is not None and point2 is not None and point3 is not None:
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point1) + np.array(point2) + np.array(point3))/4.0))

      for j in range(n_sampled_points_in_own_subdomain_y):
        #   +--3
        #  /  /*
        # +--2*1 (right center)
        # |  */
        # +--0
        for i in range(n_sampled_points_in_own_subdomain_x):
          point0 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i]
          point1 = None
          point2 = None
          point3 = None
          
          if i < n_sampled_points_in_own_subdomain_x-1:
            point1 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i+1]
          
          if k < n_sampled_points_in_own_subdomain_z-1:
            point2 = node_positions_3d_mesh[(k+1)*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i]
            
          if i < n_sampled_points_in_own_subdomain_x-1 and k < n_sampled_points_in_own_subdomain_z-1:
            point3 = node_positions_3d_mesh[(k+1)*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i+1]
            
          if point0 is not None and point2 is not None:
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point2))/2.0))
            
          if point0 is not None and point1 is not None and point2 is not None and point3 is not None:
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point1) + np.array(point2) + np.array(point3))/4.0))

        #   7--5
        #  /| *|
        # 6--*41 (center)
        # | *|/
        # 2--0
        for i in range(n_sampled_points_in_own_subdomain_x):
          point0 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i]
          point1 = None
          point2 = None
          point3 = None
          point4 = None
          point5 = None
          point6 = None
          point7 = None
          
          if i < n_sampled_points_in_own_subdomain_x-1:
            point1 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i+1]
          
          if j < n_sampled_points_in_own_subdomain_y-1:
            point2 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + (j+1)*n_sampled_points_in_own_subdomain_x + i]
            
          if i < n_sampled_points_in_own_subdomain_x-1 and j < n_sampled_points_in_own_subdomain_y-1:
            point3 = node_positions_3d_mesh[k*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + (j+1)*n_sampled_points_in_own_subdomain_x + i+1]
          
          if k < n_sampled_points_in_own_subdomain_z-1:
            point4 = node_positions_3d_mesh[(k+1)*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i]
          
          if i < n_sampled_points_in_own_subdomain_x-1 and k < n_sampled_points_in_own_subdomain_z-1:
            point5 = node_positions_3d_mesh[(k+1)*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + j*n_sampled_points_in_own_subdomain_x + i+1]
            
          if j < n_sampled_points_in_own_subdomain_y-1 and k < n_sampled_points_in_own_subdomain_z-1:
            point6 = node_positions_3d_mesh[(k+1)*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + (j+1)*n_sampled_points_in_own_subdomain_x + i]
          
          if i < n_sampled_points_in_own_subdomain_x-1 and j < n_sampled_points_in_own_subdomain_y-1 and k < n_sampled_points_in_own_subdomain_z-1:
            point7 = node_positions_3d_mesh[(k+1)*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_x + (j+1)*n_sampled_points_in_own_subdomain_x + i+1]
          
          if point0 is not None and point2 is not None and point4 is not None and point6 is not None:
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point2) + np.array(point4) + np.array(point6))/4.0))
            
          if point0 is not None and point1 is not None and point2 is not None and point3 is not None and point4 is not None and point5 is not None and point6 is not None and point7 is not None:
            node_positions_3d_mesh_quadratic.append(list((np.array(point0) + np.array(point1) + np.array(point2) + np.array(point3) + np.array(point4) + np.array(point5) + np.array(point6) + np.array(point7))/4.0))
    
    # store quadratic 3D mesh in the config
    meshes["3Dmesh_quadratic"] = {
      "nElements": variables.n_elements_3D_mesh,
      "nRanks": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z],
      "nodePositions": node_positions_3d_mesh_quadratic,
      "inputMeshIsGlobal": False,
      "setHermiteDerivatives": False,
      "logKey": "3Dmesh_quadratic"
    }
    
    # compute size for quadratic mesh, these values are only for the console output
    n_points_3D_mesh_global_x = 2*n_points_3D_mesh_global_x - 1
    n_points_3D_mesh_global_y = 2*n_points_3D_mesh_global_y - 1
    n_points_3D_mesh_global_z = 2*n_points_3D_mesh_global_z - 1
    
    n_sampled_points_in_own_subdomain_x = n_sampled_points_in_own_subdomain_x*2
    n_sampled_points_in_own_subdomain_y = n_sampled_points_in_own_subdomain_y*2
    n_sampled_points_in_own_subdomain_z = n_sampled_points_in_own_subdomain_z*2
      
    # border subdomains have one element less than variables.fibers
    if own_subdomain_coordinate_x == variables.n_subdomains_x-1:
      n_sampled_points_in_own_subdomain_x -= 1
    if own_subdomain_coordinate_y == variables.n_subdomains_y-1:
      n_sampled_points_in_own_subdomain_y -= 1
    if own_subdomain_coordinate_z == variables.n_subdomains_z-1:
      n_sampled_points_in_own_subdomain_z -= 1

  n_points_3D_mesh_global = n_points_3D_mesh_global_x*n_points_3D_mesh_global_y*n_points_3D_mesh_global_z
   
  # output for debugging
  if variables.debug_output:
    print("{}: point sampling for elements, unsampled points: {} x {} x {}, sampling stride: {}, {}, {}".format(rank_no, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber, variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z))
    print("{}: sampled points, n_points_3D_mesh_global: {} x {} x {} = sum({}) x sum({}) x sum({})".format(rank_no, n_points_3D_mesh_global_x, n_points_3D_mesh_global_y, n_points_3D_mesh_global_z, \
      [n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)],\
      [n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)],\
      [n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)]))

  # output information about partitioning on rank 0
  if rank_no == 0:
    n_states_cellml = 4
    print("{} rank{}, partitioning: x{} x y{} x z{}".format(n_ranks, "s" if n_ranks > 1 else "", variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
    print("{} x {} = {} fibers, per partition: {} x {} = {}".format(variables.n_fibers_x, variables.n_fibers_y, variables.n_fibers_total, variables.n_fibers_per_subdomain_x, variables.n_fibers_per_subdomain_y, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y))
    print("per fiber: 1D mesh nodes global: {}, local: {}".format(variables.n_points_whole_fiber, variables.n_points_per_subdomain_z))
    print(" {} 3D mesh nodes global: {} x {} x {}, local: {} x {} x {}".format("quadratic" if quadratic_3d_mesh else "   linear", n_points_3D_mesh_global_x, n_points_3D_mesh_global_y, n_points_3D_mesh_global_z, n_sampled_points_in_own_subdomain_x, n_sampled_points_in_own_subdomain_y, n_sampled_points_in_own_subdomain_z))
    print("number of degrees of freedom:")
    print("                    1D fiber: {:10d}  (per process: {})".format(variables.n_points_whole_fiber, variables.n_points_per_subdomain_z))
    print("            0D-1D monodomain: {:10d}  (per process: {})".format(variables.n_points_whole_fiber*n_states_cellml, variables.n_points_per_subdomain_z*n_states_cellml))
    print(" all fibers 0D-1D monodomain: {:10d}  (per process: {})".format(variables.n_fibers_total*variables.n_points_whole_fiber*n_states_cellml, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y*variables.n_points_per_subdomain_z*n_states_cellml))
    print("                 3D bidomain: {:10d}  (per process: {})".format(n_points_3D_mesh_global, n_sampled_points_in_own_subdomain_x*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_z))
    print("                       total: {:10d}  (per process: {})".format(variables.n_fibers_total*variables.n_points_whole_fiber*n_states_cellml+n_points_3D_mesh_global, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y*variables.n_points_per_subdomain_z*n_states_cellml+n_sampled_points_in_own_subdomain_x*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_z))

  ###############################
  # determine 1D meshes of variables.fibers

  # fiber nos of the variables.fibers that are handled on the own subdomain
  variables.fibers_on_own_rank = [fiber_no(own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y) \
    for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(own_subdomain_coordinate_y)) \
    for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(own_subdomain_coordinate_x))]
    
  if variables.debug_output:
    print("{}: rank {}, n_elements_3D_mesh: {}, subdomain coordinate ({},{},{})/({},{},{})".format(rank_no, rank_no, variables.n_elements_3D_mesh, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
    print("{}:    fibers x: [{}, {}]".format(rank_no, 0, n_fibers_in_subdomain_x(own_subdomain_coordinate_x)))
    print("{}:    fibers y: [{}, {}]".format(rank_no, 0, n_fibers_in_subdomain_y(own_subdomain_coordinate_y)))
    print("{}:       ({})".format(rank_no, variables.fibers_on_own_rank))
    print("{}:    points z: [{}, {}] ({})".format(rank_no, variables.z_point_index_start, variables.z_point_index_end, n_points_in_subdomain_z(own_subdomain_coordinate_z)))
      
  # determine number of nodes and elements of the local part of a fiber  
  variables.n_fiber_nodes_on_subdomain = n_points_in_subdomain_z(own_subdomain_coordinate_z)   # number of nodes without ghosts

  variables.fiber_start_node_no = n_points_in_previous_subdomains_z(own_subdomain_coordinate_z)

  # loop over all fibers
  for i in range(variables.n_fibers_total):

    # if fiber is computed on own rank
    if i in variables.fibers_on_own_rank:
      
      if variables.debug_output:
        print("{}: fiber {} is in fibers on own rank, {}".format(rank_no, i, str(variables.fibers_on_own_rank)))
      
      n_fiber_elements_on_subdomain = variables.n_fiber_nodes_on_subdomain

      # top subdomain has one element less than nodes
      if own_subdomain_coordinate_z == variables.n_subdomains_z-1:
        n_fiber_elements_on_subdomain -= 1

      # address fiber data
      memory_size_fiber = variables.n_points_whole_fiber * 3 * 8
      offset = 32 + header_length + i*memory_size_fiber + variables.fiber_start_node_no*3*8
      
      # fill fiber_node_positions
      # if data was loaded directly here in the python script, assign the corresponding node positions
      if load_fiber_data:
        fiber_file_handle.seek(offset)
        
        fiber_node_positions = []
        for point_no in range(variables.n_fiber_nodes_on_subdomain):
          point = []
          for j in range(3):
            double_raw = fiber_file_handle.read(8)
            value = struct.unpack('d', double_raw)[0]
            point.append(value)
          
          fiber_node_positions.append(point)
        
        if variables.debug_output:
          print("------")  
          print("i: ",i)
          print("fiber_node_positions: ",fiber_node_positions[0:10])
          print("fiber_start_node_no: ",variables.fiber_start_node_no,", n_fiber_elements_on_subdomain: ",n_fiber_elements_on_subdomain)
          print("fibers[i]: ",variables.fibers[i][variables.fiber_start_node_no:variables.fiber_start_node_no+10])
          
        if np.linalg.norm(np.array(variables.fibers[i][variables.fiber_start_node_no:variables.fiber_start_node_no+n_fiber_elements_on_subdomain]) - np.array(fiber_node_positions[0:n_fiber_elements_on_subdomain])) > 1e-3:
          print("mismatch fiber node positions!")
          quit()
          
      else:   # add command at which position the node data in the binary file can be found, the c++ core will load the data
        fiber_node_positions = [fiber_file, [(offset, variables.n_fiber_nodes_on_subdomain)]]
      
      if variables.debug_output:
        print("{}: define mesh \"{}\", n_fiber_elements_on_subdomain: {}, fiber_node_positions: {}".format(rank_no, "MeshFiber_{}".format(i), \
          str(n_fiber_elements_on_subdomain), str(fiber_node_positions)))
        
      # define 1D fiber mesh
      meshes["MeshFiber_{}".format(i)] = {
        "nElements": n_fiber_elements_on_subdomain,
        "nodePositions": fiber_node_positions,
        "inputMeshIsGlobal": False,
        "nRanks": [variables.n_subdomains_z],
        "setHermiteDerivatives": False
      }
      
    else:
      # for variables.fibers that are not computed on own rank, set empty lists for node positions and number of elements
      fiber_node_positions = []
      n_fiber_elements_on_subdomain = []

    # only add log key for fiber 0 to prevent too much data in the log files
    if i == 0 and "MeshFiber_{}".format(i) in meshes:
      meshes["MeshFiber_{}".format(i)]["logKey"] = "Fiber{}".format(i)
      
  # output more detailed partitioning information, disabled
  if rank_no == 0 and n_ranks < 10 and False:
    print("rank configuration: ")
    
    for subdomain_coordinate_y in range(variables.n_subdomains_y):
      for subdomain_coordinate_x in range(variables.n_subdomains_x):
        print("subdomain (x,y)=({},{})".format(subdomain_coordinate_x, subdomain_coordinate_y))
        print("variables.n_subdomains_z: {}".format(variables.n_subdomains_z))
        for rankNo in range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y):
          print("  rank {}".format(rankNo))
        for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
          for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
            print("  fiber {} ({},{}) in subdomain uses ranks {}".format(\
              fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y), \
              fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y, \
              list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y))))

  return [meshes, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber]
