# Multiple 1D variables.fibers (monodomain) with 3D EMG (static bidomain), biceps geometry
# arguments: -help
#
# if variables.fiber_file=cuboid.bin, it uses a small cuboid test example

import numpy as np
import pickle
import sys
import struct
import argparse
sys.path.insert(0, '..')
import variables    # file variables.py

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])    

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.own_subdomain_coordinate_x = rank_no % variables.n_subdomains_x
variables.own_subdomain_coordinate_y = (int)(rank_no / variables.n_subdomains_x) % variables.n_subdomains_y
variables.own_subdomain_coordinate_z = (int)(rank_no / variables.n_subdomains_xy)

#print("rank: {}/{}".format(rank_no,n_ranks))

# generate cuboid fiber file
if variables.fiber_file == "cuboid.bin":
  
  size_x = 1e-1
  size_y = 1e-1
  size_z = 0.2
  
  variables.n_fibers_x = 4
  variables.n_fibers_y = variables.n_fibers_x
  variables.n_points_whole_fiber = 20
  
  if rank_no == 0:
    print("create cuboid.bin with size [{},{},{}], n points [{},{},{}]".format(size_x, size_y, size_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber))
    
    # write header
    with open(variables.fiber_file, "wb") as outfile:
      
      # write header
      header_str = "opendihu self-generated cuboid  "
      outfile.write(struct.pack('32s',bytes(header_str, 'utf-8')))   # 32 bytes
      outfile.write(struct.pack('i', 40))  # header length
      outfile.write(struct.pack('i', variables.n_fibers_x*variables.n_fibers_x))   # n_fibers
      outfile.write(struct.pack('i', variables.n_points_whole_fiber))   # variables.n_points_whole_fiber
      outfile.write(struct.pack('i', 0))   # nBorderPointsXNew
      outfile.write(struct.pack('i', 0))   # nBorderPointsZNew
      outfile.write(struct.pack('i', 0))   # nFineGridFibers_
      outfile.write(struct.pack('i', 1))   # nRanks
      outfile.write(struct.pack('i', 1))   # nRanksZ
      outfile.write(struct.pack('i', 0))   # nFibersPerRank
      outfile.write(struct.pack('i', 0))   # date
    
      # loop over points
      for y in range(variables.n_fibers_x):
        for x in range(variables.n_fibers_x):
          for z in range(variables.n_points_whole_fiber):
            point = [x*(float)(size_x)/(variables.n_fibers_x-1), y*(float)(size_y)/(variables.n_fibers_y-1), z*(float)(size_z)/(variables.n_points_whole_fiber-1)]
            outfile.write(struct.pack('3d', point[0], point[1], point[2]))   # data point

# set output writer    
variables.output_writer_fibers = []
variables.output_writer_emg = []
if variables.paraview_output:
  variables.output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/emg", "binary": True, "fixedFormat": False, "combineFiles": True})
  variables.output_writer_fibers.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/variables.fibers", "binary": True, "fixedFormat": False, "combineFiles": True})

# set values for cellml model
if "shorten" in variables.cellml_file:
  variables.parameters_used_as_intermediate = [32]
  variables.parameters_used_as_constant = [65]
  variables.parameters_initial_values = [0.0, 1.0]
  variables.nodal_stimulation_current = 1200.
  
elif "hodgkin_huxley" in variables.cellml_file:
  variables.parameters_used_as_intermediate = []
  variables.parameters_used_as_constant = [2]
  variables.parameters_initial_values = [0.0]
  variables.nodal_stimulation_current = 40.

try:
  variables.fiber_file_handle = open(variables.fiber_file, "rb")
except:
  print("Error: Could not open fiber file \"{}\"".format(variables.fiber_file))
  quit()

# parse variables.fibers from a binary fiber file that was created by parallel_fiber_estimation
# parse file header to extract number of variables.fibers
bytes_raw = variables.fiber_file_handle.read(32)
header_str = struct.unpack('32s', bytes_raw)[0]
header_length_raw = variables.fiber_file_handle.read(4)
header_length = struct.unpack('i', header_length_raw)[0]

parameters = []
for i in range(int(header_length/4.) - 1):
  double_raw = variables.fiber_file_handle.read(4)
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
if variables.load_data_from_file:
  variables.fibers = []
  for fiber_no in range(variables.n_fibers_total):
    fiber = []
    for point_no in range(variables.n_points_whole_fiber):
      point = []
      for i in range(3):
        double_raw = variables.fiber_file_handle.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
      fiber.append(point)
    variables.fibers.append(fiber)
          
# callback functions
# --------------------------
def get_motor_unit_no(fiber_no):
  return int(variables.fiber_distribution[fiber_no % len(variables.fiber_distribution)]-1)

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

# for debugging output show when the first 20 variables.fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
  print("Debugging output about fiber firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(variables.firing_times,0)
  for fiber_no_index in range(variables.n_fibers_total):
    if fiber_no_index % 100 == 0:
      t_intermediate = timeit.default_timer()
      if t_intermediate - t_start > 100:
        print("Note: break after {}/{} fibers ({:.0f}%) because it already took {:.3f}s".format(fiber_no_index,variables.n_fibers_total,100.0*fiber_no_index/(variables.n_fibers_total-1.),t_intermediate - t_start))
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
  variables.fibers = []
  last_time = 0
  last_mu_no = first_stimulation_info[0][1]
  for stimulation_info in first_stimulation_info:
    mu_no = stimulation_info[1]
    fiber_no = stimulation_info[0]
    if mu_no == last_mu_no:
      variables.fibers.append(fiber_no)
    else:
      if last_time is not None:
        if len(variables.fibers) > 10:
          print("{:8.2f} {:3} {} (only showing first 10, {} total)".format(last_time,last_mu_no,str(variables.fibers[0:10]),len(variables.fibers)))
        else:
          print("{:8.2f} {:3} {}".format(last_time,last_mu_no,str(variables.fibers)))
        n_stimulated_mus += 1
      else:
        if len(variables.fibers) > 10:
          print("  never stimulated: MU {:3}, fibers {} (only showing first 10, {} total)".format(last_mu_no,str(variables.fibers[0:10]),len(variables.fibers)))
        else:
          print("  never stimulated: MU {:3}, fibers {}".format(last_mu_no,str(variables.fibers)))
        n_not_stimulated_mus += 1
      variables.fibers = [fiber_no]

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

# define helper functions for fiber numbering

# number of variables.fibers that are handled inside the subdomain x
def n_fibers_in_subdomain_x(subdomain_coordinate_x):
  a1 = variables.n_fibers_x - variables.n_subdomains_x*variables.n_fibers_per_subdomain_x              # number of subdomains with high number of variables.fibers
  a2 = variables.n_subdomains_x - a1                                               # number of subdomains with low number of variables.fibers
  if subdomain_coordinate_x < a1:
    return variables.n_fibers_per_subdomain_x+1      # high number of fibersr
  else:
    return variables.n_fibers_per_subdomain_x    # low number of variables.fibers
  
# number of variables.fibers that are handled inside the subdomain y
def n_fibers_in_subdomain_y(subdomain_coordinate_y):
  a1 = variables.n_fibers_y - variables.n_subdomains_y*variables.n_fibers_per_subdomain_y              # number of subdomains with high number of variables.fibers
  a2 = variables.n_subdomains_y - a1                                               # number of subdomains with low number of variables.fibers
  if subdomain_coordinate_y < a1:
    return variables.n_fibers_per_subdomain_y+1     # high number of variables.fibers
  else:
    return variables.n_fibers_per_subdomain_y     # low number of variables.fibers

def n_fibers_in_previous_subdomains_y(subdomain_coordinate_y):
  # number of variables.fibers handled in previous subdomains in y direction
  a1 = variables.n_fibers_y - variables.n_subdomains_y*variables.n_fibers_per_subdomain_y              # number of subdomains with high number of variables.fibers
  a2 = variables.n_subdomains_y - a1                                               # number of subdomains with low number of variables.fibers
  
  if subdomain_coordinate_y < a1:
    return subdomain_coordinate_y * (variables.n_fibers_per_subdomain_y+1)
  else:
    return a1 * (variables.n_fibers_per_subdomain_y+1) + (subdomain_coordinate_y-a1) * variables.n_fibers_per_subdomain_y
    
def n_fibers_in_previous_subdomains_x(subdomain_coordinate_x):
  # number of variables.fibers handled in previous subdomains in x direction
  a1 = variables.n_fibers_x - variables.n_subdomains_x*variables.n_fibers_per_subdomain_x              # number of subdomains with high number of variables.fibers
  a2 = variables.n_subdomains_x - a1                                               # number of subdomains with low number of variables.fibers
  
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
  a1 = variables.n_points_whole_fiber - variables.n_subdomains_z*variables.n_points_per_subdomain_z              # number of subdomains with high number of variables.fibers
  a2 = variables.n_subdomains_z - a1                                               # number of subdomains with low number of variables.fibers
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

#####################
# define 3D mesh
# determine node positions of the 3D mesh
node_positions_3d_mesh = []
if not variables.load_data_from_file:
  node_positions_3d_mesh = [variables.fiber_file, []]

# range of points in z direction
variables.z_point_index_start = n_points_in_previous_subdomains_z(variables.own_subdomain_coordinate_z)
variables.z_point_index_end = variables.z_point_index_start + n_points_in_subdomain_z(variables.own_subdomain_coordinate_z)

# loop over z point indices
for k in range(n_sampled_points_in_subdomain_z(variables.own_subdomain_coordinate_z)):
  z_point_index = variables.z_point_index_start + k*variables.sampling_stride_z
  
  if variables.own_subdomain_coordinate_z == variables.n_subdomains_z-1 and k == n_sampled_points_in_subdomain_z(variables.own_subdomain_coordinate_z)-1:
    z_point_index = variables.z_point_index_end-1
    
  #print("{}: variables.sampling_stride_z: {}, k: {}, z: {}/{}".format(rank_no, variables.sampling_stride_z, k, z_point_index, variables.z_point_index_end))
  
  # loop over variables.fibers for own rank
  # loop over fiber in y-direction
  for j in range(n_sampled_points_in_subdomain_y(variables.own_subdomain_coordinate_y)):
    fiber_in_subdomain_coordinate_y = j*variables.sampling_stride_y
    
    # on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
    if variables.own_subdomain_coordinate_y == variables.n_subdomains_y-1 and j == n_sampled_points_in_subdomain_y(variables.own_subdomain_coordinate_y)-1:
      fiber_in_subdomain_coordinate_y = n_fibers_in_subdomain_y(variables.own_subdomain_coordinate_y)-1
    
    #print("{}: variables.sampling_stride_y: {}, j: {}, y: {}/{}".format(rank_no, variables.sampling_stride_y, j, fiber_in_subdomain_coordinate_y, n_fibers_in_subdomain_y(variables.own_subdomain_coordinate_y)))
      
    # loop over fiber in x-direction
    for i in range(n_sampled_points_in_subdomain_x(variables.own_subdomain_coordinate_x)):
      fiber_in_subdomain_coordinate_x = i*variables.sampling_stride_x
      
      # on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
      if variables.own_subdomain_coordinate_x == variables.n_subdomains_x-1 and i == n_sampled_points_in_subdomain_x(variables.own_subdomain_coordinate_x)-1:
        fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(variables.own_subdomain_coordinate_x)-1
      
      #print("{}: variables.sampling_stride_x: {}, i: {}, x: {}/{}".format(rank_no, variables.sampling_stride_x, i, fiber_in_own_subdomain_coordinate_x, n_fibers_in_subdomain_x(variables.own_subdomain_coordinate_x)))
      
      # get fiber no
      fiber_index = fiber_no(variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
      
      # read point from fiber file
      memory_size_fiber = variables.n_points_whole_fiber*3*8
      offset = 32 + header_length + fiber_index*memory_size_fiber + z_point_index*3*8
      
      if variables.load_data_from_file:
        
        variables.fiber_file_handle.seek(offset)
      
        point = []
        for j in range(3):
          double_raw = variables.fiber_file_handle.read(8)
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
    n_sampled_points_in_subdomain_x(variables.own_subdomain_coordinate_x),
    n_sampled_points_in_subdomain_y(variables.own_subdomain_coordinate_y), 
    n_sampled_points_in_subdomain_z(variables.own_subdomain_coordinate_z)
  ]

# border subdomains have one element less than variables.fibers
if variables.own_subdomain_coordinate_x == variables.n_subdomains_x-1:
  variables.n_elements_3D_mesh[0] -= 1
if variables.own_subdomain_coordinate_y == variables.n_subdomains_y-1:
  variables.n_elements_3D_mesh[1] -= 1
if variables.own_subdomain_coordinate_z == variables.n_subdomains_z-1:
  variables.n_elements_3D_mesh[2] -= 1

# set the entry for the config
variables.meshes = {}
variables.meshes["3Dmesh"] = {
  "nElements": variables.n_elements_3D_mesh,
  "nRanks": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z],
  "nodePositions": node_positions_3d_mesh,
  "inputMeshIsGlobal": False,
  "setHermiteDerivatives": False,
  "logKey": "3Dmesh"
}

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
variables.linear_elasticity_dirichlet_bc = {}
for i in range(n_points_3D_mesh_global_x*n_points_3D_mesh_global_y):
  variables.linear_elasticity_dirichlet_bc[(n_points_3D_mesh_global_z-1)*n_points_3D_mesh_global_x*n_points_3D_mesh_global_y + i] = 0.0
    
# Neumann BC at bottom nodes, traction downwards
nx = n_points_3D_mesh_global_x-1
ny = n_points_3D_mesh_global_y-1
nz = n_points_3D_mesh_global_z-1
#variables.linear_elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": [0.0,0.0,-1.0], "face": "2-"} for j in range(ny) for i in range(nx)]
variables.linear_elasticity_neumann_bc = []
    
if variables.debug_output:
  print("{}: point sampling for elements, unsampled points: {} x {} x {}, sampling stride: {}, {}, {}".format(rank_no, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber, variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z))
  print("{}: sampled points, n_points_3D_mesh_global: {} x {} x {} = sum({}) x sum({}) x sum({})".format(rank_no, n_points_3D_mesh_global_x, n_points_3D_mesh_global_y, n_points_3D_mesh_global_z, \
    [n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)],\
    [n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)],\
    [n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)]))

if rank_no == 0:
  print("diffusion solver type: {}".format(variables.diffusion_solver_type))
  print("{} ranks, partitioning: x{} x y{} x z{}".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
  print("{} x {} = {} fibers, per partition: {} x {} = {}".format(variables.n_fibers_x, variables.n_fibers_y, variables.n_fibers_total, variables.n_fibers_per_subdomain_x, variables.n_fibers_per_subdomain_y, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y))
  print("{} points per fiber, per partition: {}".format(variables.n_points_whole_fiber, variables.n_points_per_subdomain_z))
  print("number of degrees of freedom:")
  print("                    1D fiber: {:10d}  (per process: {})".format(variables.n_points_whole_fiber, variables.n_points_per_subdomain_z))
  print("            0D-1D monodomain: {:10d}  (per process: {})".format(variables.n_points_whole_fiber*4, variables.n_points_per_subdomain_z*4))
  print(" all fibers 0D-1D monodomain: {:10d}  (per process: {})".format(variables.n_fibers_total*variables.n_points_whole_fiber*4, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y*variables.n_points_per_subdomain_z*4))
  print("                 3D bidomain: {:10d}  (per process: {})".format(n_points_3D_mesh_global, n_sampled_points_in_subdomain_x(variables.own_subdomain_coordinate_x)*n_sampled_points_in_subdomain_y(variables.own_subdomain_coordinate_y)*n_sampled_points_in_subdomain_z(variables.own_subdomain_coordinate_z)))
  print("                       total: {:10d}  (per process: {})".format(variables.n_fibers_total*variables.n_points_whole_fiber*4+n_points_3D_mesh_global, variables.n_fibers_per_subdomain_x*variables.n_fibers_per_subdomain_y*variables.n_points_per_subdomain_z*4+n_sampled_points_in_subdomain_x(variables.own_subdomain_coordinate_x)*n_sampled_points_in_subdomain_y(variables.own_subdomain_coordinate_y)*n_sampled_points_in_subdomain_z(variables.own_subdomain_coordinate_z)))

###############################
# determine 1D variables.meshes of variables.fibers

# fiber nos of the variables.fibers that are handled on the own subdomain
variables.fibers_on_own_rank = [fiber_no(variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y) \
  for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(variables.own_subdomain_coordinate_y)) \
  for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(variables.own_subdomain_coordinate_x))]
  
if variables.debug_output:
  print("{}: rank {}, n_elements_3D_mesh: {}, subdomain coordinate ({},{},{})/({},{},{})".format(rank_no, rank_no, variables.n_elements_3D_mesh, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z))
  print("{}:    fibers x: [{}, {}]".format(rank_no, 0, n_fibers_in_subdomain_x(variables.own_subdomain_coordinate_x)))
  print("{}:    fibers y: [{}, {}]".format(rank_no, 0, n_fibers_in_subdomain_y(variables.own_subdomain_coordinate_y)))
  print("{}:       ({})".format(rank_no, variables.fibers_on_own_rank))
  print("{}:    points z: [{}, {}] ({})".format(rank_no, variables.z_point_index_start, variables.z_point_index_end, n_points_in_subdomain_z(variables.own_subdomain_coordinate_z)))
    
# determine number of nodes and elements of the local part of a fiber  
variables.n_fiber_nodes_on_subdomain = n_points_in_subdomain_z(variables.own_subdomain_coordinate_z)   # number of nodes without ghosts

variables.fiber_start_node_no = n_points_in_previous_subdomains_z(variables.own_subdomain_coordinate_z)

# loop over all variables.fibers
for i in range(variables.n_fibers_total):

  # if fiber is computed on own rank
  if i in variables.fibers_on_own_rank:
    
    if variables.debug_output:
      print("{}: fiber {} is in fibers on own rank, {}".format(rank_no, i, str(variables.fibers_on_own_rank)))
    
    n_fiber_elements_on_subdomain = variables.n_fiber_nodes_on_subdomain

    # top subdomain has one element less than nodes
    if variables.own_subdomain_coordinate_z == variables.n_subdomains_z-1:
      n_fiber_elements_on_subdomain -= 1

    # address fiber data
    memory_size_fiber = variables.n_points_whole_fiber * 3 * 8
    offset = 32 + header_length + i*memory_size_fiber + variables.fiber_start_node_no*3*8
    
    # if data was loaded directly here in the python script, assign the corresponding node positions
    if variables.load_data_from_file:
      variables.fiber_file_handle.seek(offset)
      
      fiber_node_positions = []
      for point_no in range(variables.n_fiber_nodes_on_subdomain):
        point = []
        for j in range(3):
          double_raw = variables.fiber_file_handle.read(8)
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
          
    else:   # add command at which position the node data in the binary file can be found, the core loads the data
      fiber_node_positions = [variables.fiber_file, [(offset, variables.n_fiber_nodes_on_subdomain)]]
    
    if variables.debug_output:
      print("{}: define mesh \"{}\", n_fiber_elements_on_subdomain: {}, fiber_node_positions: {}".format(rank_no, "MeshFiber_{}".format(i), \
        str(n_fiber_elements_on_subdomain), str(fiber_node_positions)))
      
    # define mesh
    variables.meshes["MeshFiber_{}".format(i)] = {
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
  if i == 0 and "MeshFiber_{}".format(i) in variables.meshes:
    variables.meshes["MeshFiber_{}".format(i)]["logKey"] = "Fiber{}".format(i)
  
if rank_no == 0 and n_ranks < 10 and False:
  print("rank configuration: ")
  
  for subdomain_coordinate_y in range(variables.n_subdomains_y):
    for subdomain_coordinate_x in range(variables.n_subdomains_x):
      print("subdomain (x,y)=({},{})".format(subdomain_coordinate_x, subdomain_coordinate_y))
      print("n_subdomains_z: {}".format(variables.n_subdomains_z))
      for rankNo in range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y):
        print("  rank {}".format(rankNo))
      for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
          print("  fiber {} ({},{}) in subdomain uses ranks {}".format(\
            fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y), \
            fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y, \
            list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y))))

# sanity checking at the end
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
          no = fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
          if no != counter:
            print("error: fiber_no({},{},{},{}) = {}, counter = {}".format(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,no,counter))
          else:
            print("   ok: fiber_no({},{},{},{}) = {}, counter = {}".format(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,no,counter))
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
