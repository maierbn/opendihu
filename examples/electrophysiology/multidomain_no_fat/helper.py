# helper function that are used by the settings for multidomain
#
import numpy as np
import pickle
import sys
import struct
import scipy.stats

def load_mesh(fiber_file, sampling_stride_z, rank_no):

  # get the mesh nodes, either from a .bin file or a python pickle file
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
    n_points_initial_whole_fiber = parameters[1]


    # parse whole fiber file
    fiber_data = []
    mesh_node_positions = []
    for fiber_no in range(n_fibers_total):
      fiber = []
      for point_no in range(n_points_initial_whole_fiber):
        point = []
        for i in range(3):
          double_raw = fiber_file_handle.read(8)
          value = struct.unpack('d', double_raw)[0]
          point.append(value)
        fiber.append(point)
        
      # sample fiber in z direction
      new_fiber = []
      for point_no in range(0,n_points_initial_whole_fiber,sampling_stride_z):
        point = fiber[point_no]
        new_fiber.append(point)
      
      fiber_data.append(new_fiber)
              
    # set node positions
    n_points_whole_fiber = len(fiber_data[0])
    mesh_node_positions = [0 for _ in range(n_fibers_x*n_fibers_y*n_points_whole_fiber)]
    n_linear_elements_per_coordinate_direction = [n_fibers_x-1, n_fibers_y-1, n_points_whole_fiber-1]
    for k in range(n_points_whole_fiber):
      for j in range(n_fibers_y):
        for i in range(n_fibers_x):
          mesh_node_positions[k*n_fibers_x*n_fibers_y + j*n_fibers_x + i] = fiber_data[j*n_fibers_x + i][k]
          
    if rank_no == 0:
      print("fiber_file:            \"{}\"".format(fiber_file))
      print("n fibers:              {} ({} x {})".format(n_fibers_total, n_fibers_x, n_fibers_y))
      print("n points per fiber:    {}, sampling by stride {}".format(n_points_initial_whole_fiber, sampling_stride_z))
      print("3D mesh:               {} x {} x {} nodes".format(n_linear_elements_per_coordinate_direction[0]+1, n_linear_elements_per_coordinate_direction[1]+1, n_linear_elements_per_coordinate_direction[2]+1))
      
    n_points_xy = n_fibers_x*n_fibers_y
    bottom_node_indices = list(range(n_points_xy))
    n_points = n_points_xy*n_points_whole_fiber
    top_node_indices = list(range(n_points-n_points_xy,n_points))
    
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

  return (mesh_node_positions,fiber_data,bottom_node_indices,top_node_indices,n_linear_elements_per_coordinate_direction)

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
  n_fibers_x = int(np.sqrt(len(fiber_data)))

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
    point4 = np.array(fiber_data[(n_fibers_x-1)//2][z_index_fiber])
    point1 = np.array(fiber_data[n_fibers_x-1][z_index_fiber])
    point2 = np.array(fiber_data[-n_fibers_x][z_index_fiber])
    point5 = np.array(fiber_data[(-n_fibers_x)//2][z_index_fiber])
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
