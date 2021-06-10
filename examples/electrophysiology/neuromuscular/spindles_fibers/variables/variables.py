# This file contains all global variables for the fibers_emg example and their default values. These are the parameters and other internal variables.
# These values will be used by all involved scripts: helper.py, create_partitioned_meshes_for_settings.py and settings_fibers_emg.py
# settings_fibers_emg.py handles setting the parameter values. Those can be overridden on the command line and by specifying a custom variables.py script
# To run the simulation use the settings_fibers_emg.py file, which imports this file, e.g. ./fibers_emg ../settings_fibers_emg.py custom_variables.py

# material parameters
# --------------------
PMax = 7.3                          # maximum stress [N/cm^2]
Conductivity = 3.828                # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # membrane capacitance [uF/cm^2]
damping_factor = 0                  # velocity dependent damping factor

innervation_zone_width = 0.         # not used [cm], this will later be used to specify a variance of positions of the innervation point at the fibers

# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
diffusion_solver_maxit = 1e4
diffusion_solver_reltol = 1e-10
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "gamg" # preconditioner
potential_flow_solver_maxit = 1e4
potential_flow_solver_reltol = 1e-10
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False   # If the initial guess for the emg linear system should be set to the previous solution
emg_solver_maxit = 1e4
emg_solver_abstol = 1e-5
emg_solver_reltol = 1e-5

# elasticity
elasticity_solver_type = "preonly"
elasticity_preconditioner_type = "lu"
snes_max_iterations = 10                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 2       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-5      # absolute tolerance of the nonlinear solver
linear_relative_tolerance = 1e-5           # relative tolerance of the residual of the linear solver
linear_absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver

# timing parameters
# -----------------
end_time = 1000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz. This is not used here.
dt_0D = 1e-3                        # [ms] timestep width of ODEs
dt_1D = 1.5e-3                      # [ms] timestep width of diffusion
dt_multidomain = 1e-2               # [ms] timestep width of multidomain
dt_splitting = 3e-3                 # [ms] overall timestep width of strang splitting
dt_elasticity = 1e0                 # [ms] time step width of elasticity solver
output_timestep = 1e0               # [ms] timestep for output files
activation_start_time = 0           # [ms] time when to start checking for stimulation
output_timestep_multidomain = 0.1   # [ms] timestep for multidomain output files
output_timestep_elasticity = 0.1    # [ms] timestep for elasticity output files

# input files
# -----------
# CellML model, Shorten or Hodgkin-Huxley
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../../input/shorten.cpp"
#cellml_file = "../../../input/hodgkin_huxley_1952.c"

# Fiber geometry, binary file
#fiber_file = "../../../input/3000fibers.bin"
#fiber_file = "../../../input/7x7fibers.bin"
fiber_file = "../../../input/13x13fibers.bin"
#fiber_file = "../../../input/49fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"

load_fiber_data = False             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.
debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = True        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = False             # If the paraview output writer should be enabled
adios_output = False                # If the MegaMol/ADIOS output writer should be enabled
python_output = False               # If the Python output writer should be enabled
exfile_output = False               # If the Exfile output writer should be enabled
initial_guess_nonzero = True        # if the initial guess of the multidomain solver should be set to the previous values, this is only possible if an iterative solver is used
theta = 0.5                         # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True    # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False      # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
show_linear_solver_output = True    # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
optimization_type = "vc"            # the optimization_type used in the cellml adapter, "vc" uses explicit vectorization
approximate_exponential_function = False   # if the exponential function should be approximated by a Taylor series with only 11 FLOPS
dynamic = True                      # if the dynamic hyperelasticity solver should be used

# motor unit stimulation times
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
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 50
sampling_stride_fat = 1

# how much of the multidomain mesh is used for elasticity
sampling_factor_elasticity_x = 0.5    
sampling_factor_elasticity_y = 0.5
sampling_factor_elasticity_z = 0.5
sampling_factor_elasticity_fat_y = 0.5

# scenario name for log file
scenario_name = ""

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
# These functions can be redefined differently in a custom variables script
def get_am(fiber_no, mu_no):
  return Am

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  return stimulation_frequency

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return [0]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return activation_start_time


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
firing_times = None
n_fibers_per_subdomain_x = None
n_fibers_per_subdomain_y = None
n_points_per_subdomain_z = None
z_point_index_start = None
z_point_index_end = None
meshes = None
potential_flow_dirichlet_bc = None
elasticity_dirichlet_bc = None
elasticity_neumann_bc = None
fibers_on_own_rank = None
n_fiber_nodes_on_subdomain = None
fiber_start_node_no = None
generate_linear_3d_mesh = True
generate_quadratic_3d_mesh = True
fat_mesh_n_points = None
fat_mesh_n_points_global = None
local_range_i = None
local_range_k = None
relative_factors = None
n_compartments = None
rho = None
material_parameters = None
nx = None
ny = None
nz = None
constant_body_force = None
pmax = None
bottom_traction = None
cellml_file = ""
states_initial_values = []
main_elasticity_dirichlet_bc = {}
fix_bottom = False
