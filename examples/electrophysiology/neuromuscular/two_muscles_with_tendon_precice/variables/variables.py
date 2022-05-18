import os
import numpy as np

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

innervation_zone_width = 0.         # not used [cm], this will later be used to specify a variance of positions of the innervation point at the fibers

# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "none" # preconditioner
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False   # If the initial guess for the emg linear system should be set to the previous solution

# timing parameters
# -----------------
end_time = 1000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs
dt_1D = 1.5e-3                      # [ms] timestep width of diffusion
dt_splitting = 3e-3                 # [ms] overall timestep width of strang splitting
dt_3D = 1e0                         # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep = 1e0               # [ms] timestep for output files
output_timestep_3D_emg = 1e0        # [ms] timestep for output files
output_timestep_3D = 1e0            # [ms] timestep for output files
output_timestep_fibers = 1e0        # [ms] timestep for output files
activation_start_time = 0           # [ms] time when to start checking for stimulation

# input files
# -----------
input_directory = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input")

fiber_file = input_directory + "/13x13fibers.bin" #TODO: is this the muscle geometry?
fiber_distribution_file = input_directory + "/MU_fibre_distribution_3780.txt"
firing_times_file = input_directory + "/MU_firing_times_real.txt"
cellml_file = input_directory+"/hodgkin_huxley_1952.c"
muscle_spindle_cellml_file = input_directory+"/mileusenic_muscle_spindle_equilibrium.cellml"
golgi_tendon_organ_cellml_file = input_directory+"/hodgkin_huxley_1952.cellml"
motoneuron_cellml_file = input_directory+"/WSBM_1457_MN_Cisi_Kohn_2008.cellml"
interneuron_cellml_file = input_directory+"/hodgkin_huxley_1952.cellml"


load_fiber_data = False             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.
debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = True        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = False             # If the paraview output writer should be enabled
adios_output = False                # If the MegaMol/ADIOS output writer should be enabled
python_output = False               # If the Python output writer should be enabled
exfile_output = False               # If the Exfile output writer should be enabled


# motor unit stimulation times

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

mapping_tolerance = 0.1
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

muscle1_extent = [3.0, 3.0, 14.8] # [cm, cm, cm]
tendon_length = 1.2 # cm
muscle2_extent = [3.0, 3.0, 14.8] # [cm, cm, cm]

n_elements_muscle1 = [2, 2, 4] # linear elements. each qudaratic element uses the combined nodes of 8 linear elements
n_elements_muscle2 = [2, 2, 4]
n_points_whole_fiber = 40
n_fibers_x = 4
n_fibers_y = 4

# further internal variables that will be set by the helper.py script and used in the config in settings_fibers_emg.py
n_fibers_total = None
n_subdomains_xy = None
own_subdomain_coordinate_x = None
own_subdomain_coordinate_y = None
own_subdomain_coordinate_z = None
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
elasticity_dirichlet_bc = None
elasticity_neumann_bc = None
fibers_on_own_rank = None
n_fiber_nodes_on_subdomain = None
fiber_start_node_no = None
generate_linear_3d_mesh = False
generate_quadratic_3d_mesh = True
rho = None
material_parameters = None
nx = None
ny = None
nz = None
constant_body_force = None
pmax = None
bottom_traction = None
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1
states_initial_values = []
enable_coupling = True
enable_force_length_relation = True
lambda_dot_scaling_factor = 1


# muscle spindles
n_muscle_spindles = 3
muscle_spindle_mappings = {
  # first: don't define any mappings. then fix all warnings and errors

  # define the order of the parameters for initial values and callbacks
  # opendihu                    cellml
  ("parameter", 0):            "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("parameter", 1):            "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("parameter", 2):            "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("parameter", 3):            "modell/gamma_sta",                 # Static fusimotor drive
  ("parameter", 4):            "modell/gamma_dyn",                 # Dynamic fusimotor drive
  # we first have to define the cellml constants as parameters. primary_afferent ist schon ein state. 
  # opendihu                    cellml
  ("connectorSlot", "ms_out"): "modell/primary_afferent",          # Ia Afferent firing frequency
  ("connectorSlot", "ms_in0"): "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("connectorSlot", "ms_in1"): "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("connectorSlot", "ms_in2"): "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("connectorSlot", "ms_in3"): "modell/gamma_sta",                 # Static fusimotor drive
  ("connectorSlot", "ms_in4"): "modell/gamma_dyn",                 # Dynamic fusimotor drive
}
muscle_spindle_parameters_initial_values = [0, 0, 0, 5, 0]    # [L, L_dot, L_ddot, gamma_sta, gamma_dyn] (see above)
muscle_spindle_delay = 30             # [ms] signal delay between muscle spindle model and motoneuron model


# golgi tendon organs
n_golgi_tendon_organs = 3
golgi_tendon_organ_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("connectorSlot", "gt_out"): "membrane/V",        # voltage
  ("connectorSlot", "gt_in"):  "membrane/i_Stim",   # stimulation
}
golgi_tendon_organ_parameters_initial_values = [0]    # [i_Stim]
golgi_tendon_organ_delay = 300


# load cortical input values
cortical_input_file = input_directory+"/cortical_input_realistic.txt"
cortical_input = np.genfromtxt(cortical_input_file, delimiter=",")


# motor neurons
n_motoneurons = 101
n_motor_units = n_motoneurons
motoneuron_mappings = {
  ("parameter", 0):            "motor_neuron/drive",   # stimulation
  ("parameter", 1):            "lumped_geometry_parameters/C_m",
  ("parameter", 2):            "lumped_geometry_parameters/R_i",
  ("parameter", 3):            "lumped_geometry_parameters/R_md",
  ("parameter", 4):            "lumped_geometry_parameters/R_ms",
  ("parameter", 5):            "lumped_geometry_parameters/l_d",
  ("parameter", 6):            "lumped_geometry_parameters/l_s",
  ("parameter", 7):            "lumped_geometry_parameters/r_d",
  ("parameter", 8):            "lumped_geometry_parameters/r_s",
  ("connectorSlot", "mn_out"): "motor_neuron/V_s",     # voltage
  ("connectorSlot", "mn_in"):  "motor_neuron/drive",   # stimulation
}
# initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
motoneuron_parameters_initial_values = []
# Parameter ranges for MN Paramters [smallest MN value, largest MN value]
C_m = [1, 1]
Ri = [0.07, 0.07]
r_d = [20.75e-4, 46.25e-4]
l_d = [0.55, 1.06]
r_s = [38.75e-4, 56.5e-4]
l_s = [77.5e-4, 113e-4]
R_ms = [1.15, 0.65]
R_md = [14.4, 6.05]
gg = np.linspace(0,1,num = n_motoneurons)
#
for mu_idx in range(n_motoneurons):
  # compute a scaling factor that distributes variables exponentially between a lower and an upper value
  factor = np.exp(np.log(100)*gg[mu_idx-1])/100  # 100^(gg-1) in [0,1]
  # add parameter values for motoneuron mu_idx
  motoneuron_parameters_initial_values += [0.0, C_m[0]+factor*(C_m[1]-C_m[0]), Ri[0]+factor*(Ri[1]-Ri[0]), 
                                          R_md[0]+factor*(R_md[1]-R_md[0]), R_ms[0]+factor*(R_ms[1]-R_ms[0]), 
                                          l_d[0]+factor*(l_d[1]-l_d[0]), l_s[0]+factor*(l_s[1]-l_s[0]), 
                                          r_d[0]+factor*(r_d[1]-r_d[0]), r_s[0]+factor*(r_s[1]-r_s[0])]

# interneurons
n_interneurons = 3
interneuron_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("parameter", 1):            "membrane/Cm",       # stimulation
  ("connectorSlot", "in_out"): "membrane/V",        # voltage
  ("connectorSlot", "in_in"):  "membrane/i_Stim",   # stimulation
}
# initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
interneuron_parameters_initial_values = []    # [i_Stim, Cm]
for i in range(n_interneurons):
  # compute a scaling factor that runs exponentially from min_factor to max_factor
  min_factor = 0.5
  max_factor = 2.0
  
  # ansatz scaling_factor(i) = c1 + c2*exp(i),
  # scaling_factor(0) = min = c1 + c2  =>  c1 = min - c2
  # scaling_factor(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_factor - min_factor) / (np.exp(n_interneurons-1) - 1)
  c1 = min_factor - c2
  scaling_factor = c1 + c2*np.exp(i)
  
  # add parameter values for motoneuron i
  interneuron_parameters_initial_values += [0.0, 1*scaling_factor]

