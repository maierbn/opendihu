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
        mesh_node_positions.append(point)
      
      fiber_data.append(new_fiber)
              
    # set node positions
    n_points_whole_fiber = len(fiber_data[0])
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

  return (mesh_node_positions,fiber_data,bottom_node_indices,top_node_indices,n_linear_elements_per_coordinate_direction)


def compute_compartment_relative_factors(mesh_node_positions, fiber_data, motor_units):
    
  # list of fibers, fiber = list of points, point = list with 3 coordinate entries
  n_compartments = len(motor_units)

  # create relative factors for compartments
  #if rank_no == 0:
  #  print("determine relative factors for {} motor units:\n{}".format(n_compartments, motor_units))

  # create data structure with 0
  relative_factors = np.zeros((n_compartments, len(mesh_node_positions)))   # each row is one compartment

  # loop over nodes of mesh
  for node_no,node_position in enumerate(mesh_node_positions):
    node_position = np.array(node_position)
    
    # loop over motor units
    for motor_unit_no,motor_unit in enumerate(motor_units):
      
      # find point on fiber that is closest to current node
      fiber_no = motor_unit["fiber_no"]
      if fiber_no >= len(fiber_data):
        print("Error with motor unit {}, only {} fibers available".format(motor_unit, len(fiber_data)))
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
  
  return relative_factors
