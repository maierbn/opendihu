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
  a1 = (int)((variables.n_fibers_x - variables.n_subdomains_x*variables.n_fibers_per_subdomain_x) / variables.granularity_x)     # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_x - a1                                      # number of subdomains with low number of fibers
  if subdomain_coordinate_x < a1:
    return variables.n_fibers_per_subdomain_x + variables.granularity_x   # high number of fibers
  elif subdomain_coordinate_x < variables.n_subdomains_x-1:
    return variables.n_fibers_per_subdomain_x                             # low number of fibers
  else:
    return variables.n_fibers_per_subdomain_x + variables.n_fibers_x % variables.granularity_x   # last subdomain has low number of fibers + granularity remainder
    
# number of fibers that are handled inside the subdomain y
def n_fibers_in_subdomain_y(subdomain_coordinate_y):
  a1 = (int)((variables.n_fibers_y - variables.n_subdomains_y*variables.n_fibers_per_subdomain_y) / variables.granularity_y)    # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_y - a1                                      # number of subdomains with low number of fibers
  if subdomain_coordinate_y < a1:
    return variables.n_fibers_per_subdomain_y + variables.granularity_y   # high number of fibers
  elif subdomain_coordinate_y < variables.n_subdomains_y-1:
    return variables.n_fibers_per_subdomain_y                             # low number of fibers
  else:
    return variables.n_fibers_per_subdomain_y + variables.n_fibers_y % variables.granularity_y  # last subdomain has low number of fibers + 1 (only for granularity_y==2)

def n_fibers_in_previous_subdomains_x(subdomain_coordinate_x):
  # number of fibers handled in previous subdomains in x direction
  a1 = (int)((variables.n_fibers_x - variables.n_subdomains_x*variables.n_fibers_per_subdomain_x) / variables.granularity_x)     # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_x - a1                                               # number of subdomains with low number of fibers
  
  if subdomain_coordinate_x < a1:
    return subdomain_coordinate_x * (variables.n_fibers_per_subdomain_x+variables.granularity_x)
  else:
    return a1 * (variables.n_fibers_per_subdomain_x+variables.granularity_x) + (subdomain_coordinate_x-a1) * variables.n_fibers_per_subdomain_x

def n_fibers_in_previous_subdomains_y(subdomain_coordinate_y):
  # number of fibers handled in previous subdomains in y direction
  a1 = (int)((variables.n_fibers_y - variables.n_subdomains_y*variables.n_fibers_per_subdomain_y) / variables.granularity_y)    # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_y - a1                                               # number of subdomains with low number of fibers
  
  if subdomain_coordinate_y < a1:
    return subdomain_coordinate_y * (variables.n_fibers_per_subdomain_y+variables.granularity_y)
  else:
    return a1 * (variables.n_fibers_per_subdomain_y+variables.granularity_y) + (subdomain_coordinate_y-a1) * variables.n_fibers_per_subdomain_y
    
# global fiber no, from subdomain coordinate and coordinate inside the subdomain
def get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y):
  return (n_fibers_in_previous_subdomains_y(subdomain_coordinate_y) + fiber_in_subdomain_coordinate_y)*variables.n_fibers_x \
    + n_fibers_in_previous_subdomains_x(subdomain_coordinate_x) + fiber_in_subdomain_coordinate_x

# number of points that are handled inside the subdomain z (without ghost points)
def n_points_in_subdomain_z(subdomain_coordinate_z):
  a1 = (int)((variables.n_points_whole_fiber - variables.n_subdomains_z*variables.n_points_per_subdomain_z) / variables.granularity_z)  # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_z - a1                                      # number of subdomains with low number of fibers
  if subdomain_coordinate_z < a1:
    return variables.n_points_per_subdomain_z + variables.granularity_z     # high number of points
  elif subdomain_coordinate_z < variables.n_subdomains_z-1:
    return variables.n_points_per_subdomain_z                             # low number of points
  else:
    return variables.n_points_per_subdomain_z + variables.n_points_whole_fiber % variables.granularity_z # last subdomain has low number of fibers + 1
  
def n_points_in_previous_subdomains_z(subdomain_coordinate_z):
  # number of points handled in previous subdomains in z direction
  a1 = (int)((variables.n_points_whole_fiber - variables.n_subdomains_z*variables.n_points_per_subdomain_z) / variables.granularity_z)  # number of subdomains with high number of fibers
  a2 = variables.n_subdomains_z - a1                                               # number of subdomains with low number of points
  
  if subdomain_coordinate_z < a1:
    return subdomain_coordinate_z * (variables.n_points_per_subdomain_z+variables.granularity_z)
  else:
    return a1 * (variables.n_points_per_subdomain_z+variables.granularity_z) + (subdomain_coordinate_z-a1) * variables.n_points_per_subdomain_z

# number of points in 3D mesh in the subdomain, in x direction
def n_sampled_points_in_subdomain_x(subdomain_coordinate_x):
  n = n_fibers_in_subdomain_x(subdomain_coordinate_x)       # there are as many "linear" elements as nodes
  if subdomain_coordinate_x == variables.n_subdomains_x-1:  # only at the right boundary there is one more point, i.e. one element less than points
    n -= 1
  
  if variables.generate_quadratic_3d_mesh:        # if a quadratic mesh is requested
    result = (int)(np.floor(n / (variables.sampling_stride_x*2)) * 2)   # make sure the number of elements in the subdomain is even for quadratic elements (this has nothing to do with granularity)
  else:                                           # if a linear mesh is requested
    result = (int)(np.ceil(n / variables.sampling_stride_x))            # for linear elements, make remainder elements smaller (->ceil), for quadratic elements, make them larger (->floor)
    
  if subdomain_coordinate_x == variables.n_subdomains_x-1:
    result += 1
  return result

# number of points in 3D mesh in the subdomain, in y direction
def n_sampled_points_in_subdomain_y(subdomain_coordinate_y):
  n = n_fibers_in_subdomain_y(subdomain_coordinate_y)       # there are as many "linear" elements as nodes
  if subdomain_coordinate_y == variables.n_subdomains_y-1:  # only at the right boundary there is one more point, i.e. one element less than points
    n -= 1
  
  if variables.generate_quadratic_3d_mesh:        # if a quadratic mesh is requested
    result = (int)(np.floor(n / (variables.sampling_stride_y*2)) * 2)   # make sure the number of elements in the subdomain is even for quadratic elements (this has nothing to do with granularity)
  else:                                           # if a linear mesh is requested
    result = (int)(np.ceil(n / variables.sampling_stride_y))            # for linear elements, make remainder elements smaller (->ceil), for quadratic elements, make them larger (->floor)
    
  if subdomain_coordinate_y == variables.n_subdomains_y-1:
    result += 1
  return result

# number of points in 3D mesh in the subdomain, in z direction
def n_sampled_points_in_subdomain_z(subdomain_coordinate_z):
  n = n_points_in_subdomain_z(subdomain_coordinate_z)       # there are as many "linear" elements as nodes
  if subdomain_coordinate_z == variables.n_subdomains_z-1:  # only at the right boundary there is one more point, i.e. one element less than points
    n -= 1
  
  if variables.generate_quadratic_3d_mesh:        # if a quadratic mesh is requested
    result = (int)(np.floor(n / (variables.sampling_stride_z*2)) * 2)   # make sure the number of elements in the subdomain is even for quadratic elements (this has nothing to do with granularity)
  else:                                           # if a linear mesh is requested
    result = (int)(np.ceil(n / variables.sampling_stride_z))            # for linear elements, make remainder elements smaller (->ceil), for quadratic elements, make them larger (->floor)
    
  if subdomain_coordinate_z == variables.n_subdomains_z-1:
    result += 1
  return result

# this is the main function that creates the meshes and partitioning
def create_partitioned_meshes_for_settings(n_subdomains_x, n_subdomains_y, n_subdomains_z, 
                                           fiber_file, load_fiber_data, sampling_stride_x, sampling_stride_y, sampling_stride_z, 
                                           generate_linear_3d_mesh=False, generate_quadratic_3d_mesh=False, fiber_set_rank_nos=False, 
                                           have_fibers=True, include_global_node_positions=False):
  """
  Parse the binary fiber geometry file and creates 1D fiber meshes and a 3D mesh (linear or quadratic elements for the 3D mesh)
  It create the `meshes` Dict that can directly be used in the Python config like: 
  config = {
    "Meshes": meshes,
    ...
  }
  
  One or two 3D meshes with keys "3Dmesh" and "3Dmesh_quadratic" will be created, and multiple 1D fiber meshes with keys "MeshFiber_0", "MeshFiber_1", etc.
  The 3D meshes are aligned with the 1D meshes.
  
  The mesh specifications have the following form:
  A 1D mesh:
  meshes["MeshFiber_0"] = {
    "nElements":             n_fiber_elements_on_subdomain,    # the number of 1D elements of this fiber on the own processes' subdomain
    "nodePositions":         fiber_node_positions,             # all node positions of this fiber on the own processes' subdomain
    "inputMeshIsGlobal":     False,                            # `False` means that only locale nodes and elements are provided
    "nRanks":                [variables.n_subdomains_z],       # the number of ranks in z-direction over which this fiber mesh is spread
    "setHermiteDerivatives": False,
    #"rankNos":              [0,1],                            # the ranks of this fiber, only set if argument fiber_set_rank_nos is true
  }
  
  The 3D mesh (if generate_linear_3d_mesh):
  meshes["3Dmesh"] = {
    "nElements":             variables.n_elements_3D_mesh_linear,     # the number of elements on the own process
    "nRanks":                [n_subdomains_x, n_subdomains_y, n_subdomains_z],   # the number of subdomains = processes of the 3D mesh
    "nodePositions":         node_positions_3d_mesh,           # all local node positions
    "inputMeshIsGlobal":     False,                            # `False` means that only locale nodes and elements are provided
    "setHermiteDerivatives": False,
    "logKey":                "3Dmesh",
    "globalNodePositions":   [...],                            # only if option include_global_node_positions is True
  }
  and/or (if generate_quadratic_3d_mesh):
  meshes["3Dmesh_quadratic"] = {
    "nElements":             variables.n_elements_3D_mesh_quadratic,             # the number of elements on the own process
    "nRanks":                [n_subdomains_x, n_subdomains_y, n_subdomains_z],   # the number of subdomains = processes of the 3D mesh
    "nodePositions":         node_positions_3d_mesh_quadratic,                   # all local node positions
    "inputMeshIsGlobal":     False,                                              # `False` means that only locale nodes and elements are provided
    "setHermiteDerivatives": False,
    "logKey": "3Dmesh_quadratic",
    
    # information on how many nodes there are in the quadratic 3D mesh, this is not needed for the opendihu core 
    "nPointsLocal": [n_sampled_points_in_own_subdomain_x, n_sampled_points_in_own_subdomain_y, n_sampled_points_in_own_subdomain_z],
    "nPointsGlobal": [n_points_3D_mesh_global_x, n_points_3D_mesh_global_y, n_points_3D_mesh_global_z],
  }
  
  The quadratic 3D mesh: same nodes as linear 3D mesh, but under meshes["3Dmesh_quadratic"], number of elements is smaller by factor of 8, 
  but has quadratic elements with 27 nodes instead of 8.
  
  This means that all quantities like, e.g. the local 3D node positions can be retrieved by meshes["3Dmesh"]["nodePositions"] and so on.
  
  :param n_subdomains_x:  number of subdomains in x direction.
  :param n_subdomains_y:  number of subdomains in y direction
  :param n_subdomains_z:  number of subdomains in z direction, i.e. the number of subdivisions per fiber. 
                          The z direction is the direction of the fibers along the muscle.
                          the product n_subdomains_x*n_subdomains_y*n_subdomains_z has to equal the number of processes
  :param fiber_file: file name of the fiber file "*.bin" which contains fiber geometry in 3D space. The file has a pickle format.
  :param load_fiber_data: If the actual geometry data should be read from the file and the node positions placed in the generated meshes Dict. 
                          If this is true, the node positions will be in the mesh config.
                          If this is false, it inserts the filename and the position/offset in the file where the geometry data is present.
                          This information is later used by the C++ code to read the file in parallel. This is necessary for very large runs,
                          where the file has to be parsed in parallel and not already here in the python script.
  :param sampling_stride_x:   Grid point stride in x direction to make the 3D mesh coarser than the grid points of the 1D fibers.
                              E.g. 2 means there will be 2 fibers per 3D element in x-direction on average.
                              (Or 3 if you want, two of them are on the edges of the 3D elements, one is along the center)
  :param sampling_stride_y:   Grid point stride in y direction to make the 3D mesh coarser than the grid points of the 1D fibers.
  :param sampling_stride_z:   Grid point stride in z direction to make the 3D mesh coarser than the grid points of the 1D fibers.
                              E.g. 20 means that there will be 20 mesh nodes of the 1D mesh per 3D element. 
                              The z direction is the direction of the fibers along the muscle.
  :param generate_linear_3d_mesh: Whether to create the linear mesh (under "3Dmesh")
  :param generate_quadratic_3d_mesh: Whether to create the quadratic mesh (under "3Dmesh_quadratic")
  :param fiber_set_rank_nos:  If the "rankNos" option of the fiber meshes should be set to the ranks that take part in the computation of the fiber.
  :param have_fibers:         If fiber meshes should be created.
  :param include_global_node_positions: If the global node positions for the "3Dmesh" should be added under the key "globalNodePositions"
  
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
  variables.generate_quadratic_3d_mesh = generate_quadratic_3d_mesh
  variables.generate_linear_3d_mesh = generate_linear_3d_mesh
  
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
    print("\033[0;31mError: Could not open fiber file \"{}\".\033[0m".format(fiber_file))
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
  
  # The number of points in z direction in the input file is `variables.n_points_whole_fiber`.
  variables.n_points_whole_fiber = parameters[1]
  
  # if the old option local_sampling_stride_z instead of sampling_stride_z is set, print a warning
  if hasattr(variables, 'local_sampling_stride_z'):
    if variables.sampling_stride_z == 1 and variables.local_sampling_stride_z != 1:
      variables.sampling_stride_z = variables.local_sampling_stride_z
      if rank_no == 0:
        print("Note, local_sampling_stride_z is no longer used! Simply use sampling_stride_z! (Now setting sampling_stride_z = {})".format(variables.sampling_stride_z))
  
  if "cuboid.bin" not in fiber_file:  
    variables.n_fibers_x = (int)(np.round(np.sqrt(variables.n_fibers_total)))
    variables.n_fibers_y = variables.n_fibers_x

  if rank_no == 0 and have_fibers:
    print("n fibers:              {} ({} x {}), sampled by stride {} x {}".format(variables.n_fibers_total, variables.n_fibers_x, variables.n_fibers_y, variables.sampling_stride_x, variables.sampling_stride_y))
    print("n points per fiber:    {}, sampled by stride {}".format(variables.n_points_whole_fiber, variables.sampling_stride_z))
    
  # granularity_x,granularity_y,granularity_z are the lowest granularity of the number of fibers per subdomain
  # setting these variables to a higher value than 1 leads to subdomains that have to include a number of elements that is a multiple of granularity 
  variables.granularity_x = 1   
  variables.granularity_y = 1
  variables.granularity_z = 1
    
  if not hasattr(variables, 'distribute_nodes_equally'):
    variables.distribute_nodes_equally = False
    
  if not variables.distribute_nodes_equally:
    variables.granularity_x = variables.sampling_stride_x
    variables.granularity_y = variables.sampling_stride_y
    variables.granularity_z = variables.sampling_stride_z
  
  # if only one of linear or quadratic mesh should be created, it is created directly
  # if both the quadratic and the linear mesh should be created, the quadratic mesh is created first and then the linear mesh is created from the quadratic
    
  # if a quadratic mesh is required
  if generate_quadratic_3d_mesh:
    #load_fiber_data = True
    variables.granularity_x = 2
    variables.granularity_y = 2
    variables.granularity_z = 2
    
    if not variables.distribute_nodes_equally:
      variables.granularity_x = max(2, variables.sampling_stride_x // 2 * 2)    # the granularity for quadratic meshes has to be a multiple of 2 and at least 2
      variables.granularity_y = max(2, variables.sampling_stride_y // 2 * 2)
      variables.granularity_z = max(2, variables.sampling_stride_z // 2 * 2)
    
    if variables.n_fibers_x % 2 == 0:
      print("\033[0;31mError: Quadratic mesh is requested but number of fibers in x direction ({}) is even. It has to be odd.\033[0m".format(variables.n_fibers_x))
      quit()
    if variables.n_fibers_y % 2 == 0:
      print("\033[0;31mError: Quadratic mesh is requested but number of fibers in y direction ({}) is even. It has to be odd.\033[0m".format(variables.n_fibers_y))
      quit()
    if variables.n_points_whole_fiber % 2 == 0:
      print("\033[0;31mError: Quadratic mesh is requested but number of points for fiber in z direction ({}) is even. It has to be odd.\033[0m".format(variables.n_points_whole_fiber, variables.sampling_stride_z))
      quit()
      
  # parse whole fiber file, only if enabled
  if load_fiber_data or include_global_node_positions:
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
      print("\n\n\033[1;31;40mError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n\033[0m".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
      quit()
    
  # compute average number of points per subdomain, the actual number is equal to this or +(variables.granularity-1)
  variables.n_fibers_per_subdomain_x = (int)((variables.n_fibers_x / variables.n_subdomains_x) // variables.granularity_x * variables.granularity_x)
  variables.n_fibers_per_subdomain_y = (int)((variables.n_fibers_y / variables.n_subdomains_y) // variables.granularity_y * variables.granularity_y)
  variables.n_points_per_subdomain_z = (int)((variables.n_points_whole_fiber / variables.n_subdomains_z) // variables.granularity_z * variables.granularity_z)

  if variables.n_fibers_per_subdomain_x == 0 or variables.n_fibers_per_subdomain_y == 0:
    print("\033[0;31mError: Cannot partition {}x{} fibers into {}x{} subdomains.\033[0m".format(variables.n_fibers_x, variables.n_fibers_y, variables.n_subdomains_x, variables.n_subdomains_y))
    quit()
  if variables.n_points_per_subdomain_z == 0:
    print("\033[0;31mError: Cannot partition {} points per fibers into {} subdomains.\033[0m".format(variables.n_points_per_subdomain_z, variables.n_subdomains_z))
    quit()

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
    
    if own_subdomain_coordinate_z == variables.n_subdomains_z-1:
      if k == n_sampled_points_in_own_subdomain_z-1:
        z_point_index = variables.z_point_index_end-1
      elif k == n_sampled_points_in_own_subdomain_z-2 and k >= 1:
        # for quadratic meshes, set second last point at the center between third last point and last point
        z_point_index = int(0.5*(variables.z_point_index_start + (k-1)*variables.sampling_stride_z + variables.z_point_index_end-1))
        
      
    #print("{}: sampling_stride_z: {}, k: {}/{}, z: {}/{}".format(rank_no, variables.sampling_stride_z, k, n_sampled_points_in_own_subdomain_z, z_point_index, variables.z_point_index_end))
    
    # loop over fibers for own rank
    # loop over fiber in y-direction
    for j in range(n_sampled_points_in_own_subdomain_y):
      fiber_in_subdomain_coordinate_y = j*variables.sampling_stride_y
      
      # on boundary rank set last node positions to be the boundary nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
      if own_subdomain_coordinate_y == variables.n_subdomains_y-1:
        if j == n_sampled_points_in_own_subdomain_y-1:
          fiber_in_subdomain_coordinate_y = n_fibers_in_subdomain_y(own_subdomain_coordinate_y)-1
        elif j == n_sampled_points_in_own_subdomain_y-2 and j >= 1:
          # for quadratic meshes, set second last point at the center between third last point and last point
          fiber_in_subdomain_coordinate_y = int(0.5*((j-1)*variables.sampling_stride_y + n_fibers_in_subdomain_y(own_subdomain_coordinate_y)-1))
      
      #if k==0:
      #  print("{}: sampling_stride_y: {}, j: {}/{}, y: {}/{}".format(rank_no, variables.sampling_stride_y, j, n_sampled_points_in_own_subdomain_y, fiber_in_subdomain_coordinate_y, n_fibers_in_subdomain_y(own_subdomain_coordinate_y)))
        
      # loop over fiber in x-direction
      for i in range(n_sampled_points_in_own_subdomain_x):
        fiber_in_subdomain_coordinate_x = i*variables.sampling_stride_x
        
        # on boundary rank set last node positions to be the boundary nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
        if own_subdomain_coordinate_x == variables.n_subdomains_x-1:
          if i == n_sampled_points_in_own_subdomain_x-1:
            fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(own_subdomain_coordinate_x)-1
          elif i == n_sampled_points_in_own_subdomain_x-2 and i >= 1:
            # for quadratic meshes, set second last point at the center between third last point and last point
            fiber_in_subdomain_coordinate_x = int(0.5*((i-1)*variables.sampling_stride_x + n_fibers_in_subdomain_x(own_subdomain_coordinate_x)-1))
        
        #if j == 0 and k == 0:
        #  print("{}: sampling_stride_x: {}, i: {}/{}, x: {}/{} (j={},k={})".format(rank_no, variables.sampling_stride_x, i, n_sampled_points_in_own_subdomain_x, fiber_in_subdomain_coordinate_x, n_fibers_in_subdomain_x(own_subdomain_coordinate_x),j,k))
        
        # get fiber no
        fiber_index = get_fiber_no(own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
        
        # read point from fiber file
        memory_size_fiber = variables.n_points_whole_fiber*3*8
        offset = 32 + header_length + fiber_index*memory_size_fiber + z_point_index*3*8
        
        if load_fiber_data:
          
          fiber_file_handle.seek(offset)
        
          point = []
          for component_no in range(3):
            double_raw = fiber_file_handle.read(8)
            value = struct.unpack('d', double_raw)[0]
            point.append(value)
        
          if fiber_index >= len(variables.fibers):
            print("Error: fiber_index: {}, n fibers loaded: {}, subdomain_coordinate: ({},{}), fiber coordinate in subdomain: ({},{})".format(fiber_index, len(variables.fibers),own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y))
          if z_point_index >= len(variables.fibers[fiber_index]):
            print("Error: z_point_index: {}, n points in fiber: {}".format(z_point_index, len(variables.fibers[fiber_index])))
        
          reference_point = variables.fibers[fiber_index][z_point_index]
          
          difference = np.linalg.norm(np.array(reference_point) - np.array(point))
          if difference > 1e-3:
            print("\033[0;31mError, point does not match: reference_point: ", reference_point, ", point: ", point, "\033[0m")
            quit()
          node_positions_3d_mesh.append(point)
            
          #if j == n_sampled_points_in_own_subdomain_y-1 and k == 0:
          #  print("{}: sampling_stride_x: {}, i: {}, x: {}/{} (j={},k={}) point: {} (f{}, z{})".format(rank_no, variables.sampling_stride_x, i, fiber_in_subdomain_coordinate_x, n_fibers_in_subdomain_x(own_subdomain_coordinate_x),j,k,point,fiber_index,z_point_index))
          
        else:
          node_positions_3d_mesh[1].append((offset, 1))   # command to read 1 point from offset
        #print("{}: i={},j={},k={}, point: {}".format(rank_no, i, j,k,point))
       
       
  # set local number of elements for the 3D mesh
  variables.n_elements_3D_mesh_linear = [
      n_sampled_points_in_own_subdomain_x,
      n_sampled_points_in_own_subdomain_y, 
      n_sampled_points_in_own_subdomain_z
    ]
  
  # boundary subdomains have one element less than fibers
  if own_subdomain_coordinate_x == variables.n_subdomains_x-1:
    variables.n_elements_3D_mesh_linear[0] -= 1
  if own_subdomain_coordinate_y == variables.n_subdomains_y-1:
    variables.n_elements_3D_mesh_linear[1] -= 1
  if own_subdomain_coordinate_z == variables.n_subdomains_z-1:
    variables.n_elements_3D_mesh_linear[2] -= 1

  variables.n_points_3D_mesh_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
  variables.n_points_3D_mesh_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
  variables.n_points_3D_mesh_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
  
  # set the entry for the config
  meshes = {}
  if generate_linear_3d_mesh:
    meshes["3Dmesh"] = {
      "nElements": variables.n_elements_3D_mesh_linear,
      "nRanks": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z],
      "nodePositions": node_positions_3d_mesh,      # without ghosts
      "inputMeshIsGlobal": False,
      "setHermiteDerivatives": False,
      "logKey": "3Dmesh",
      
      # set information on how many nodes there are in the quadratic 3D mesh, this is not needed for the opendihu core but might be useful in some settings script
      "nPointsLocal": [n_sampled_points_in_own_subdomain_x, n_sampled_points_in_own_subdomain_y, n_sampled_points_in_own_subdomain_z],
      "nPointsGlobal": [variables.n_points_3D_mesh_global_x, variables.n_points_3D_mesh_global_y, variables.n_points_3D_mesh_global_z],
    }
    
    # add all global node positions of all ranks
    if include_global_node_positions:
      
      meshes["3Dmesh"]["globalNodePositions"] = []
      
      # ---- loop over subdomains and nodes in subdomains in z direction ----
      for subdomain_coordinate_z in range(variables.n_subdomains_z):
        z_point_index_start = n_points_in_previous_subdomains_z(subdomain_coordinate_z)
        z_point_index_end = z_point_index_start + n_points_in_subdomain_z(subdomain_coordinate_z)
        
        # loop over z point indices
        for k in range(n_sampled_points_in_subdomain_z(subdomain_coordinate_z)):
          z_point_index = z_point_index_start + k*variables.sampling_stride_z
          if subdomain_coordinate_z == variables.n_subdomains_z-1:
            if k == n_sampled_points_in_subdomain_z(subdomain_coordinate_z)-1:
              z_point_index = z_point_index_end-1
            elif k == n_sampled_points_in_subdomain_z(subdomain_coordinate_z)-2 and k >= 1:
              # for quadratic meshes, set second last point at the center between third last point and last point
              z_point_index = int(0.5*(z_point_index_start + (k-1)*variables.sampling_stride_z + z_point_index_end-1))
          
          # ---- loop over subdomains and nodes in subdomains in y direction ----
          for subdomain_coordinate_y in range(variables.n_subdomains_y):
            # loop over fiber in y-direction
            for j in range(n_sampled_points_in_subdomain_y(subdomain_coordinate_y)):
              fiber_in_subdomain_coordinate_y = j*variables.sampling_stride_y
              # on boundary subdomain set last node positions to be the boundary nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
              if subdomain_coordinate_y == variables.n_subdomains_y-1:
                if j == n_sampled_points_in_subdomain_y(subdomain_coordinate_y)-1:
                  fiber_in_subdomain_coordinate_y = n_fibers_in_subdomain_y(subdomain_coordinate_y)-1
                elif j == n_sampled_points_in_subdomain_y(subdomain_coordinate_y)-2 and j >= 1:
                  # for quadratic meshes, set second last point at the center between third last point and last point
                  fiber_in_subdomain_coordinate_y = int(0.5*((j-1)*variables.sampling_stride_y + n_fibers_in_subdomain_y(subdomain_coordinate_y)-1))
              
              # ---- loop over subdomains and nodes in subdomains in x direction ----
              for subdomain_coordinate_x in range(variables.n_subdomains_x):
                # loop over fiber in x-direction
                for i in range(n_sampled_points_in_subdomain_x(subdomain_coordinate_x)):
                  fiber_in_subdomain_coordinate_x = i*variables.sampling_stride_x
                  # on boundary subdomain set last node positions to be the boundary nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
                  if subdomain_coordinate_x == variables.n_subdomains_x-1:
                    if i == n_sampled_points_in_subdomain_x(subdomain_coordinate_x)-1:
                      fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(subdomain_coordinate_x)-1
                    elif i == n_sampled_points_in_subdomain_x(subdomain_coordinate_x)-2 and i >= 1:
                      # for quadratic meshes, set second last point at the center between third last point and last point
                      fiber_in_subdomain_coordinate_x = int(0.5*((i-1)*variables.sampling_stride_x + n_fibers_in_subdomain_x(subdomain_coordinate_x)-1))
                      
                  # get fiber no
                  fiber_index = get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
                  
                  point = variables.fibers[fiber_index][z_point_index]
            
                  meshes["3Dmesh"]["globalNodePositions"].append(point)
  
  # create quadratic 3D mesh
  if generate_quadratic_3d_mesh:
    
    if variables.n_elements_3D_mesh_linear[0] % 2 != 0 or variables.n_elements_3D_mesh_linear[1] % 2 != 0 or variables.n_elements_3D_mesh_linear[2] % 2 != 0:
      print("\033[0;31mError, local number of elements (n_elements_3D_mesh_linear) is not even as needed for quadratic elements: {}\033[0m".format(variables.n_elements_3D_mesh_linear))
      quit()
    
    variables.n_elements_3D_mesh_quadratic = [(int)(variables.n_elements_3D_mesh_linear[0]/2), (int)(variables.n_elements_3D_mesh_linear[1]/2), (int)(variables.n_elements_3D_mesh_linear[2]/2)]
    
    # store quadratic 3D mesh in the config
    meshes["3Dmesh_quadratic"] = {
      "nElements": variables.n_elements_3D_mesh_quadratic,
      "nRanks": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z],
      "nodePositions": node_positions_3d_mesh,
      "inputMeshIsGlobal": False,
      "setHermiteDerivatives": False,
      "logKey": "3Dmesh_quadratic",
      
      # set information on how many nodes there are in the quadratic 3D mesh, this is not needed for the opendihu core but might be useful in some settings script
      "nPointsLocal": [n_sampled_points_in_own_subdomain_x, n_sampled_points_in_own_subdomain_y, n_sampled_points_in_own_subdomain_z],
      "nPointsGlobal": [variables.n_points_3D_mesh_global_x, variables.n_points_3D_mesh_global_y, variables.n_points_3D_mesh_global_z],
    }
    
  n_points_3D_mesh_global = variables.n_points_3D_mesh_global_x*variables.n_points_3D_mesh_global_y*variables.n_points_3D_mesh_global_z
   
  # output for debugging
  if variables.debug_output:
    print("{}: point sampling for elements, unsampled points: {} x {} x {}, sampling stride: {}, {}, {}".format(rank_no, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber, variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z))
    print("{}: sampled points, n_points_3D_mesh_global: {} x {} x {} = sum({}) x sum({}) x sum({})".format(rank_no, variables.n_points_3D_mesh_global_x, variables.n_points_3D_mesh_global_y, variables.n_points_3D_mesh_global_z, \
      [n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)],\
      [n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)],\
      [n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)]))
      
    print("{}: own subdomain coordinates: ({},{},{})/({},{},{})".format(rank_no, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
    print("{}: 3Dmesh           nElements: {} nRanks: {} len(nodePositions): {} ".format(rank_no, variables.n_elements_3D_mesh_linear, meshes["3Dmesh"]["nRanks"], len(node_positions_3d_mesh)))
    if "3Dmesh_quadratic" in meshes:
      print("{}: 3Dmesh_quadratic nElements: {} nRanks: {}".format(rank_no, meshes["3Dmesh_quadratic"]["nElements"], meshes["3Dmesh_quadratic"]["nRanks"]))
    

  # compute helper variables for output and checking if partitioning is valid
  if generate_quadratic_3d_mesh:
    n_elements_3D_quadratic_global_x = (int)((variables.n_points_3D_mesh_global_x-1)/2)
    n_elements_3D_quadratic_global_y = (int)((variables.n_points_3D_mesh_global_y-1)/2)
    n_elements_3D_quadratic_global_z = (int)((variables.n_points_3D_mesh_global_z-1)/2)
    n_elements_3D_quadratic_local = variables.n_elements_3D_mesh_quadratic[0] * variables.n_elements_3D_mesh_quadratic[1] * variables.n_elements_3D_mesh_quadratic[2]
    n_elements_3D_quadratic_global = n_elements_3D_quadratic_global_x * n_elements_3D_quadratic_global_y * n_elements_3D_quadratic_global_z
  if generate_linear_3d_mesh:
    n_elements_3D_linear_global_x = (int)(variables.n_points_3D_mesh_global_x-1)
    n_elements_3D_linear_global_y = (int)(variables.n_points_3D_mesh_global_y-1)
    n_elements_3D_linear_global_z = (int)(variables.n_points_3D_mesh_global_z-1)
    n_elements_3D_linear_local = variables.n_elements_3D_mesh_linear[0] * variables.n_elements_3D_mesh_linear[1] * variables.n_elements_3D_mesh_linear[2]
    n_elements_3D_linear_global = n_elements_3D_linear_global_x * n_elements_3D_linear_global_y * n_elements_3D_linear_global_z
  variables.n_subdomains = variables.n_subdomains_x * variables.n_subdomains_y * variables.n_subdomains_z
      
  # output information about partitioning on rank 0
  if rank_no == 0:      
    n_states_cellml = 4
    if "shorten" in variables.cellml_file:
      n_states_cellml = 57
    elif "slow_TK_2014" in variables.cellml_file:
      n_states_cellml = 56
    elif "hodgkin_huxley-razumova" in variables.cellml_file:
      n_states_cellml = 9

    print("{} rank{}, partitioning: x{} x y{} x z{}".format(n_ranks, "s" if n_ranks > 1 else "", variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
    if have_fibers:
      if n_ranks > 0:
        print("{} x {} = {} fibers, per partition: {} x {} = {}".format(variables.n_fibers_x, variables.n_fibers_y, variables.n_fibers_total, variables.n_fibers_per_subdomain_x, variables.n_fibers_per_subdomain_y, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y))
      else:
        print("{} x {} = {} fibers".format(variables.n_fibers_x, variables.n_fibers_y, variables.n_fibers_total))
      print("per fiber: 1D mesh    nodes global: {}, local: {}".format(variables.n_points_whole_fiber, n_points_in_subdomain_z(own_subdomain_coordinate_z)))
    
    print("  sampling 3D mesh with stride {} x {} x {} {}".format(variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z, "\n  distribute_nodes_equally: True" if variables.distribute_nodes_equally else ""))
    if generate_linear_3d_mesh:
      print("    linear 3D mesh    nodes global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
        variables.n_points_3D_mesh_global_x, variables.n_points_3D_mesh_global_y, variables.n_points_3D_mesh_global_z, n_points_3D_mesh_global, 
        n_sampled_points_in_own_subdomain_x, n_sampled_points_in_own_subdomain_y, n_sampled_points_in_own_subdomain_z, n_sampled_points_in_own_subdomain_x*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_z))
      print("    linear 3D mesh elements global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
        n_elements_3D_linear_global_x, n_elements_3D_linear_global_y, n_elements_3D_linear_global_z, n_elements_3D_linear_global,
        variables.n_elements_3D_mesh_linear[0], variables.n_elements_3D_mesh_linear[1], variables.n_elements_3D_mesh_linear[2], n_elements_3D_linear_local))
    if generate_quadratic_3d_mesh:
      print(" quadratic 3D mesh    nodes global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
        variables.n_points_3D_mesh_global_x, variables.n_points_3D_mesh_global_y, variables.n_points_3D_mesh_global_z, n_points_3D_mesh_global, 
        n_sampled_points_in_own_subdomain_x, n_sampled_points_in_own_subdomain_y, n_sampled_points_in_own_subdomain_z, n_sampled_points_in_own_subdomain_x*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_z))
      print(" quadratic 3D mesh elements global: {} x {} x {} = {}, local: {} x {} x {} = {}".format(
        n_elements_3D_quadratic_global_x, n_elements_3D_quadratic_global_y, n_elements_3D_quadratic_global_z, n_elements_3D_quadratic_global,
        variables.n_elements_3D_mesh_quadratic[0], variables.n_elements_3D_mesh_quadratic[1], variables.n_elements_3D_mesh_quadratic[2], n_elements_3D_quadratic_local))
        
    if have_fibers:
      print("number of degrees of freedom:")
      print("                    1D fiber: {:10d}  (per process: {})".format(variables.n_points_whole_fiber, n_points_in_subdomain_z(own_subdomain_coordinate_z)))
      print("            0D-1D monodomain: {:10d}  (per process: {})".format(variables.n_points_whole_fiber*n_states_cellml, n_points_in_subdomain_z(own_subdomain_coordinate_z)*n_states_cellml))
      print(" all fibers 0D-1D monodomain: {:10d}  (per process: {})".format(variables.n_fibers_total*variables.n_points_whole_fiber*n_states_cellml, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y*n_points_in_subdomain_z(own_subdomain_coordinate_z)*n_states_cellml))
      print("                 3D bidomain: {:10d}  (per process: {})".format(n_points_3D_mesh_global, n_sampled_points_in_own_subdomain_x*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_z))
      print("                       total: {:10d}  (per process: {})".format(variables.n_fibers_total*variables.n_points_whole_fiber*n_states_cellml+n_points_3D_mesh_global, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y*n_points_in_subdomain_z(own_subdomain_coordinate_z)*n_states_cellml+n_sampled_points_in_own_subdomain_x*n_sampled_points_in_own_subdomain_y*n_sampled_points_in_own_subdomain_z))

  # exit if number of elements is <= 0 on any rank
  if generate_quadratic_3d_mesh:
    if variables.n_elements_3D_mesh_quadratic[0] <= 0 or variables.n_elements_3D_mesh_quadratic[1] <= 0 or variables.n_elements_3D_mesh_quadratic[2] <= 0:
      print("\n\033[0;31mError! When partitioning {}x{}x{} quadratic 3D elements to {}x{}x{}={} ranks, rank {} gets {}x{}x{}={} elements (subdomain coordinates (0-based): ({},{},{})/({},{},{})).\n (Sampling 3D mesh with stride {} x {} x {})\nDecrease number of processes or increase mesh size or specify different partitioning.\n\033[0m".
      format(n_elements_3D_quadratic_global_x, n_elements_3D_quadratic_global_y, n_elements_3D_quadratic_global_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains,
      rank_no, variables.n_elements_3D_mesh_quadratic[0], variables.n_elements_3D_mesh_quadratic[1], variables.n_elements_3D_mesh_quadratic[2], n_elements_3D_quadratic_local,
      own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z,
      variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z))
      quit()
  if generate_linear_3d_mesh:
    if variables.n_elements_3D_mesh_linear[0] <= 0 or variables.n_elements_3D_mesh_linear[1] <= 0 or variables.n_elements_3D_mesh_linear[2] <= 0:
      print("\n\033[0;31mError! When partitioning {}x{}x{} 3D elements to {}x{}x{}={} ranks, rank {} gets {}x{}x{}={} elements (subdomain coordinates (0-based): ({},{},{})/({},{},{})).\n (Sampling 3D mesh with stride {} x {} x {})\nDecrease number of processes or increase mesh size or specify different partitioning.\n\033[0m".
      format(n_elements_3D_linear_global_x, n_elements_3D_linear_global_y, n_elements_3D_linear_global_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains,
      rank_no, variables.n_elements_3D_mesh_linear[0], variables.n_elements_3D_mesh_linear[1], variables.n_elements_3D_mesh_linear[2], n_elements_3D_linear_local,
      own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z,
      variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z))
      quit()
    
  ###############################
  # determine 1D meshes of fibers

  if have_fibers:

    # fiber nos of the fibers that are handled on the own subdomain
    variables.fibers_on_own_rank = [get_fiber_no(own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y) \
      for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(own_subdomain_coordinate_y)) \
      for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(own_subdomain_coordinate_x))]
      
    if variables.debug_output:
      print("{}: rank {}, n_elements_3D_mesh_linear: {}, subdomain coordinate ({},{},{})/({},{},{})".format(rank_no, rank_no, variables.n_elements_3D_mesh_linear, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
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
      if i in variables.fibers_on_own_rank or fiber_set_rank_nos:
        
        if variables.debug_output:
          if variables.fibers_on_own_rank:
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
          
          if variables.debug_output and False:
            print("------")  
            print("i: ",i)
            print("fiber_node_positions: ",fiber_node_positions[0:10])
            print("fiber_start_node_no: ",variables.fiber_start_node_no,", n_fiber_elements_on_subdomain: ",n_fiber_elements_on_subdomain)
            print("fibers[i]: ",variables.fibers[i][variables.fiber_start_node_no:variables.fiber_start_node_no+10])
            
          if np.linalg.norm(np.array(variables.fibers[i][variables.fiber_start_node_no:variables.fiber_start_node_no+n_fiber_elements_on_subdomain]) - np.array(fiber_node_positions[0:n_fiber_elements_on_subdomain])) > 1e-3:
            print("\033[0;31mmismatch fiber node positions!\033[0m")
            quit()
            
        else:   # add command at which position the node data in the binary file can be found, the c++ core will load the data
          fiber_node_positions = [fiber_file, [(offset, variables.n_fiber_nodes_on_subdomain)]]
        
        if variables.debug_output and False:
          print("{}: define mesh \"{}\", n_fiber_elements_on_subdomain: {}, fiber_node_positions: {}".format(rank_no, "MeshFiber_{}".format(i), \
            str(n_fiber_elements_on_subdomain), str(fiber_node_positions)))
          
        # define 1D fiber mesh
        meshes["MeshFiber_{}".format(i)] = {
          "nElements": n_fiber_elements_on_subdomain,
          "nodePositions": fiber_node_positions,
          "inputMeshIsGlobal": False,
          "nRanks": [variables.n_subdomains_z],
          "setHermiteDerivatives": False,
        }
        
        if fiber_set_rank_nos:
          rank_nos = []
          fiber_coordinate_x = i % variables.n_fibers_x
          fiber_coordinate_y = (int)(i / variables.n_fibers_y)
          
          for k in range(variables.n_subdomains_z):
            
            for x in range(variables.n_subdomains_x):
              if n_fibers_in_previous_subdomains_x(x) + n_fibers_in_subdomain_x(x) > fiber_coordinate_x:
                break
            
            for y in range(variables.n_subdomains_y):
              if n_fibers_in_previous_subdomains_y(y) + n_fibers_in_subdomain_y(y) > fiber_coordinate_y:
                break
            
            rank_no_fiber = k*variables.n_subdomains_x*variables.n_subdomains_y + y*variables.n_subdomains_x + x
            rank_nos.append(rank_no_fiber)
          
          meshes["MeshFiber_{}".format(i)]["rankNos"] = rank_nos
          print("{}: fiber {} gets rank_nos {}".format(rank_no,i,rank_nos))
        
      else:
        # for fibers that are not computed on own rank, set empty lists for node positions and number of elements
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
                get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y), \
                fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y, \
                list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y))))

  return [meshes, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber]
