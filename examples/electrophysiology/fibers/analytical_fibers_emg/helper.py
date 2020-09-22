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

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    variables.fiber_file, variables.load_fiber_data,
    variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z, variables.generate_linear_3d_mesh, variables.generate_quadratic_3d_mesh)
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result
  
variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# set output writer    
variables.output_writer_emg = []

subfolder = ""
if variables.paraview_output:
  if variables.adios_output:
    subfolder = "paraview/"
  variables.output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_big), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
  
if variables.adios_output:
  if variables.paraview_output:
    subfolder = "adios/"
  variables.output_writer_emg.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "useFrontBackBuffer": False, "combineNInstances": 1, "fileNumbering": "incremental"})
  
if variables.python_output:
  if variables.adios_output:
    subfolder = "python/"
  variables.output_writer_emg.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "binary": True, "fileNumbering": "incremental"})
  
if variables.exfile_output:
  if variables.adios_output:
    subfolder = "exfile/"
  variables.output_writer_emg.append({"format": "Exfile", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg", "fileNumbering": "incremental"})
  
# callback functions
# --------------------------
def get_motor_unit_no(fiber_no):
  return int(variables.fiber_distribution[fiber_no % len(variables.fiber_distribution)]-1)

def get_diffusion_prefactor(fiber_no, mu_no):
  return variables.get_conductivity(fiber_no, mu_no) / (variables.get_am(fiber_no, mu_no) * variables.get_cm(fiber_no, mu_no))

def last_stimulation_time_of_fiber(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(fiber_no)*alpha)
  
  # determine last firing time of fiber
  index = int(np.ceil(current_time * frequency))
  n_firing_times = np.size(variables.firing_times,0)
  
  current_index = index % n_firing_times
  firing_indices = sorted(np.where(variables.firing_times[0:current_index, mu_no] == 1))


  last_firing_index = None
  if np.size(firing_indices) > 0:
    last_firing_index = firing_indices[0][-1]
    
  if not last_firing_index:
  #  print("t:{}, i:{}, firing_times: {}, firing_indices: {}, last_firing_index: {}".format(current_time, current_index, variables.firing_times[0:current_index, mu_no], firing_indices, last_firing_index))
    return None
    
  last_firing_time = last_firing_index / frequency
  #print("t:{}, i:{}, firing_times: {}, firing_indices: {}, last_firing_index: {}, last_firing_time: {}".format(current_time, current_index, variables.firing_times[0:current_index, mu_no], firing_indices, last_firing_index, last_firing_time))
  
  return last_firing_time

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
  
# callback function for transmembrane voltage, Vm, instead of monodomain
def set_vm_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, fiber_no):
  # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
  # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
  #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
  # time_step_no:        (int)   current time step number
  # current_time:        (float) the current simulation time
  # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
  #                       i.e. [point0_component0, point0_component1, ... pointN_component0, point1_component0, point1_component1, ...]
  #                       After the call, these values will be assigned to the field variable.
  # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
  # additional_argument: The value of the option "additionalArgument", can be any Python object.
  
  # This function computes the Vm values for a (subdomain of a) single fiber with given fiber_no
  
  def g(z):
    """ Rosenfalck phenomenological model """
    if z >= 0:
      return 96 * z**3 * np.exp(-z) - 90
    else:
      return -90
  
  # determine when the current fiber was stimulated the last time, according to motor unit and firing_times_file
  last_stimulation_time = last_stimulation_time_of_fiber(fiber_no, variables.stimulation_frequency, current_time)

  # loop over nodes of local partition in fiber
  for local_dof_no in range(len(values)):
    
    # get the global no. of the current local dof no.
    global_dof_no = global_natural_dofs[local_dof_no]
    
    # compute distance k from center N/2 of fiber
    n_nodes_per_fiber = n_nodes_global_per_coordinate_direction[0]
    N = n_nodes_per_fiber
    k = global_dof_no
    if k > N/2:
      k = k - N/2
    else:
      k = N/2 - k
    
    # transform position on fiber from cm to mm
    k = k / variables.n_nodes_per_mm      # current position on the fiber in [mm]
    propagation_velocity = 4    # parameter conduction velocity [mm/ms]
    
    # if there was no stimulation yet, use equilibrium value -90
    if last_stimulation_time is None:
      values[local_dof_no] = -90
    else:
      # if there was a stimulation, at last_stimulation_time, use Rosenfalck model
      values[local_dof_no] = g(-k + propagation_velocity*(current_time - last_stimulation_time))
      
      # debugging output
      #  print("  g(-{} + v*({}-{})) = g({}) = {}".format(k, current_time, last_stimulation_time, -k + propagation_velocity*(current_time - last_stimulation_time), values[local_dof_no]))

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
  
