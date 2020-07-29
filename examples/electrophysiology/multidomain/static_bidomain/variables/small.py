
# scenario name for log file, also subdirectory for output files
scenario_name = "small"

# material parameters
# --------------------

# timing and activation parameters
# -----------------

end_time = 40.0                    # [ms] end time of the simulation
dt_coupling = 1e-1                     # [ms] timestep for coupling between static bidomain solver and prescribed values
output_timestep_3D_emg = 1          # [ms] timestep for output files for 3D EMG

# input files
fiber_file = "../../../input/left_biceps_brachii_3x3fibers.bin"
#fiber_file = "../../../input/left_biceps_brachii_9x9fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 50

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = False

import numpy as np

# callback function for artifical activation values, instead of monodomain
def set_artifical_activation_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, additional_argument):
  # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
  # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
  #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
  # time_step_no:        (int)   current time step number
  # current_time:        (float) the current simulation time
  # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
  #                       i.e. [point0_component0, point0_component1, ... pointN_component0, point1_component0, point1_component1, ...]
  # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
  # additional_argument: The value of the option "additionalArgument", can be any Python object.

  # number of nodes in x, y and z direction
  mx = n_nodes_global_per_coordinate_direction[0]
  my = n_nodes_global_per_coordinate_direction[1]
  mz = n_nodes_global_per_coordinate_direction[2]
  
  # loop over local dofs
  for local_dof_no in range(len(values)):
    global_dof_no = global_natural_dofs[local_dof_no]
    
    # get position in the mesh of the current node
    i = global_dof_no % mx                    # index in x direction
    j = int((global_dof_no % (mx*my)) / mx)   # index in y direction
    k = int(global_dof_no / (mx*my))          # index in z direction
    
    # set activation value ∈ [0,1]
    t = current_time
    T = 10.0
    values[local_dof_no] = np.sin((0.25*float(k)/mz+t/T) * 2*np.pi) ** 2

