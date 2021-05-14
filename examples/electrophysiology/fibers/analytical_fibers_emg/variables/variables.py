# This file contains all global variables for the fibers_emg example and their default values. These are the parameters and other internal variables.
# These values will be used by all involved scripts: helper.py, create_partitioned_meshes_for_settings.py and settings_fibers_emg.py
# settings_fibers_emg.py handles setting the parameter values. Those can be overridden on the command line and by specifying a custom variables.py script
# To run the simulation use the settings_fibers_emg.py file, which imports this file, e.g. ./fibers_emg ../settings_fibers_emg.py custom_variables.py

# solvers
# -------
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "none" # preconditioner
potential_flow_solver_maxit = 1e4
potential_flow_solver_reltol = 1e-10
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False   # If the initial guess for the emg linear system should be set to the previous solution
emg_solver_maxit = 1e4
emg_solver_reltol = 1e-5

# timing parameters
# -----------------
end_time = 1000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_3D = 1e0                         # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep = 1e0               # [ms] timestep for output files
output_timestep_fibers = 1e0        # [ms] timestep for output files of fibers
output_timestep_big = 1e0           # [ms] timestep for output files of 3D intramuscular EMG data

# input files
# -----------
# Fiber geometry, binary file
fiber_file = "../../../input/left_biceps_brachii_13x13fibers.bin"

load_fiber_data = False             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.
debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = False        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = False             # If the paraview output writer should be enabled
adios_output = False                # If the MegaMol/ADIOS output writer should be enabled
python_output = False               # If the Python output writer should be enabled
exfile_output = False               # If the Exfile output writer should be enabled
enable_surface_emg = False          # Enables the surface emg output writer
enable_matrix_output = False        # if the matrix should be output to a MATLAB file

# motor unit stimulation times
cellml_file = ""
fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../../input/MU_firing_times_real.txt"
#firing_times_file = "../../../input/MU_firing_times_immediately.txt"

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 50

n_nodes_per_mm = 10

# scenario name for log file
scenario_name = ""

# further internal variables that will be set by the helper.py script and used in the config in settings_fibers_emg.py
n_fibers_total = None
n_subdomains_xy = None
own_subdomain_coordinate_x = None
own_subdomain_coordinate_y = None
own_subdomain_coordinate_z = None
n_fibers_x = None
n_fibers_y = None
n_points_whole_fiber = None
n_points_3D_mesh_global_x = None
n_points_3D_mesh_global_y = None
n_points_3D_mesh_global_z = None
output_writer_fibers = None
output_writer_emg = None
output_writer_0D_states = None
states_output = False
parameters_used_as_algebraic = None
parameters_used_as_constant = None
parameters_initial_values = None
output_algebraic_index = None
output_state_index = None
nodal_stimulation_current = None
fiber_file_handle = None
fibers = None
fiber_distribution = None
firing_times = None
n_fibers_per_subdomain_x = None
n_fibers_per_subdomain_y = None
n_points_per_subdomain_z = None
z_point_index_start = None
z_point_index_end = None
meshes = None
potential_flow_dirichlet_bc = None
use_elasticity_dirichlet_bc = None
use_elasticity_neumann_bc = None
fibers_on_own_rank = None
n_fiber_nodes_on_subdomain = None
fiber_start_node_no = None
generate_linear_3d_mesh = True
generate_quadratic_3d_mesh = False
mappings = None
enable_weak_scaling = False
custom_meshes = {}
