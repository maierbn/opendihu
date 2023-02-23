# This file contains all global variables for the fibers_emg example and their default values. These are the parameters and other internal variables.
# These values will be used by all involved scripts: helper.py, create_partitioned_meshes_for_settings.py and settings_fibers_emg.py
# settings_fibers_emg.py handles setting the parameter values. Those can be overridden on the command line and by specifying a custom variables.py script
# To run the simulation use the settings_fibers_emg.py file, which imports this file, e.g. ./fibers_emg ../settings_fibers_emg.py custom_variables.py

# material parameters
# --------------------
Pmax = 7.3                          # maximum stress [N/cm^2]
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
end_time = 20.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz. This is not used here.
dt_neuron_system            = 1e-3  # [ms]
dt_muscle_spindles          = 1e-3  # [ms]
dt_golgi_tendon_organs      = 1e-3  # [ms]
dt_interneuron              = 1e-3  # [ms]
dt_motoneuron               = 1e-3  # [ms]
dt_0D = 0.5e-3                        # [ms] timestep width of ODEs
dt_1D = 1e-3                      # [ms] timestep width of diffusion
dt_bidomain = 1e-2                  # [ms] timestep width of multidomain
dt_splitting_0D1D = 1e-3            # [ms] overall timestep width of strang splitting
dt_elasticity = 1e0                 # [ms] time step width of elasticity solver
output_timestep = 1e0               # [ms] timestep for output files
activation_start_time = 0           # [ms] time when to start checking for stimulation
output_timestep_fibers = 2   # [ms] timestep for multidomain output files
output_timestep_elasticity = 1    # [ms] timestep for elasticity output files
output_timestep_emg = 20    # [ms] timestep for emg output files

output_timestep_golgi_tendon_organs = 20
output_timestep_spindles = 1         # [ms] timestep for output of files for all sensor organs and neurons
output_timestep_motoneuron = 1    # [ms] timestep for output of files for motoneuron
output_timestep_interneurons = 1  # [ms] timestep for output of files for intermotoneuron
output_timestep_surface = 20



# input files
# -----------
# CellML model, Shorten or Hodgkin-Huxley
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../../input/shorten.cpp"
#cellml_file = "../../../input/hodgkin_huxley_1952.c"

debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = True        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = True             # If the paraview output writer should be enabled
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
#firing_times_file = "../../../input/MU_firing_times_real_no_firing.txt" # no firing

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


muscle1_extent = [3.0, 3.0, 14.8] # [cm, cm, cm]
n_elements_muscle1 = [2, 2, 4] # linear elements. each qudaratic element uses the combined nodes of 8 linear elements

n_points_whole_fiber = 40
n_fibers_x = 4
n_fibers_y = 4



## currently undefined 
maximum_number_of_threads = 1
use_aovs_memory_layout = True
enable_surface_emg = False
fast_monodomain_solver_optimizations = True


# REMOVE ? 
fiber_file = None # create_partitioned_meshes_for_settings not needed? reads fibers
load_fiber_data = False
include_global_node_positions = False

# further internal variables that will be set by the helper.py script and used in the config in settings_fibers_emg.py
n_fibers_total = None
n_subdomains_xy = None
own_subdomain_coordinate_x = 0 # TODO fix this for parallelization
own_subdomain_coordinate_y = 0 # TODO fix this for parallelization
own_subdomain_coordinate_z = 0 # TODO fix this for parallelization
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
meshes = {}
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
nx = None
ny = None
nz = None
constant_body_force = None
bottom_traction = None
states_initial_values = []
fix_bottom = False




# general parameters
# -----------------------------
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)
# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2
c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter
# for debugging, b = 0 leads to normal Mooney-Rivlin

material_parameters = [c1, c2, b, d]   # material parameters


#--------------------------------

import os
import sys
import numpy as np

input_directory = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input")

#cellml_file = input_directory+"/hodgkin_huxley_1952.c"
cellml_file = input_directory+"/hodgkin_huxley-razumova.cellml"
fiber_distribution_file = input_directory+"/MU_fibre_distribution_multidomain_67x67_100.txt"
 

# neurons and sensors
# -------------------

# muscle spindles
n_muscle_spindles = 3
muscle_spindle_cellml_file = input_directory+"/mileusenic_muscle_spindle_equilibrium.cellml"
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
golgi_tendon_organ_cellml_file = input_directory+"/hodgkin_huxley_1952.cellml"
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
motoneuron_cellml_file = input_directory+"/WSBM_1457_MN_Cisi_Kohn_2008.cellml"
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
interneuron_cellml_file = input_directory+"/hodgkin_huxley_1952.cellml"
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



def get_from_obj(data, path):
    for elem in path:
        if type(elem) == str:
            data = data[elem]
        elif type(elem) == int:
            data = data[elem]
        elif type(elem) == tuple:
            # search for key == value with (key, value) = elem
            key, value = elem
            data = next(filter(lambda e: e[key] == value, data))
        else:
            raise KeyError(f"Unknown type of '{elem}': '{type(elem)}'. Path: '{'.'.join(path)}'")
    return data




def callback_muscle_spindles_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # map from λ in the 3D mesh to muscle spindles model input
  # input_values contains a list of λ values from the muscle spindle nodes
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_muscle_spindles (per output slot if there are multiple)
  
  print()
  print()
  print("callback_muscle_spindles_input, input={}".format(input_values), flush=True)
  
  # initialize buffer, buffer is needed to compute velocity and acceleration
  if "stretch" not in buffer:
    buffer["stretch"]      = [1 for _ in range(n_input_values)]
    buffer["velocity"]     = [0 for _ in range(n_input_values)]
    buffer["acceleration"] = [0 for _ in range(n_input_values)]
    buffer["current_time"] = 0
    delta_t = 1
  else:  
    delta_t  = buffer["current_time"] - current_time
    
  # normal stimulation: parse input value of λ
  for i in range(n_input_values):
    stretch = input_values[i]
    
    # the first time when the stretch is not yet computed it has a value of 0, set to 1
    if stretch == 0:
      stretch = 1
      
    # compute velocity and acceleration by difference quotient
    velocity     = (buffer["stretch"][i] - stretch) / delta_t
    acceleration = (buffer["velocity"][i] - velocity) / delta_t
    buffer["stretch"][i]      = stretch
    buffer["velocity"][i]     = velocity
    buffer["acceleration"][i] = acceleration
    print("spindle no. {:2d}, stretch: {:+.8e}, v: {:+.8e}, a: {:+.8e}".format(i, stretch, velocity, acceleration))
      
    # scale value for muscle spindle model 
    output_values[0][i] = stretch                 # fiber stretch, output units: [L0] with L0=cm
    output_values[1][i] = velocity*1e-3           # velocity, output units: [cm*s^-1], 1e-3 cm*ms^-1 = 1e-3 cm*(1e-3s)^-1 = 1e-3 cm*1e3*s^-1 = 1 cm*s^-1
    output_values[2][i] = acceleration*1e-6       # acceleration, output units: [cm*s^-2], 1e-6 cm*ms^-2 = 1e-6 cm*(1e-3s)^-2 = 1e-6 cm*(1e6)*s^-2 = 1 cm*s^-2
    output_values[3][i] = 0                       # static fusimotor drive
    output_values[4][i] = 0                       # dynamic fusimotor drive
  
  # store current time
  buffer["current_time"] = current_time

  sys.stdout.flush() # flush output. neccessary if stdout ist a pipe (eg | tee). files seem to be ok
  
  # artifical muscle spindle input for debugging
  if False:
    for i in range(n_output_values):
      T = 1000   # [ms] cycle duration of stimulus 
      output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 0.001
      output_values[0][i] = 1
      
      # "ms_in0" -> "modell/L",                         # Spinndle Stretch (= fiber stretch)
      # "ms_in1" -> "modell/L_dot",                     # Spinndle velocity (= dL/dt)
      # "ms_in2" -> "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
      # "ms_in3" -> "modell/gamma_sta",                 # Static fusimotor drive
      # "ms_in4" -> "modell/gamma_dyn",                 # Dynamic fusimotor drive
      output_values[0][i] = 1 + np.sin(current_time / T * np.pi) * 0.1                  # Stretch    
      output_values[1][i] = np.cos(current_time / T * np.pi) * 0.1 * (np.pi/T)          # Velocity
      output_values[2][i] = np.sin(current_time / T * np.pi) * (-0.1) * (np.pi/T)**2    # Acceleration
      output_values[3][i] = 0     # static fusimotor drive
      output_values[4][i] = 0     # dynamic fusimotor drive


def callback_muscle_spindles_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping muscle spindles output -> motor neuron signals, delay signals by muscle_spindle_delay
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_motoneurons

  # initialize buffer the first time, buffer later stores time of last activation, to model signal delay
  if 0 not in buffer:
    for muscle_spindle_index in range(n_input_values):
      buffer[muscle_spindle_index] = None
  
  # total summed up signal
  total_signal = 0
      
  # loop over muscle spindles
  for muscle_spindle_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[muscle_spindle_index] > 0:
      buffer[muscle_spindle_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[muscle_spindle_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = muscle_spindle_delay             # [ms] delay of the signal
      gaussian_std_dev = 10                      # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[muscle_spindle_index]) * 5
        
      # sum up all input signals
      total_signal += delayed_signal * 1e-3
    
      print(" spindle {:2d}/{:2d} signal: {:.2e}".format(muscle_spindle_index,n_input_values,delayed_signal * 1e-3))
        
  # compute average signal over all muscle spindles
  total_signal /= n_input_values
    
  # motor neuron fires with ~14Hz if drive(t) = 5e-3
  
  # set values to all connected motoneurons
  for motoneuron_index in range(n_output_values):
    output_values[0][motoneuron_index] = total_signal
  
  print("muscle_spindles_to_motoneurons: {} -> {}".format(input_values, output_values))
  sys.stdout.flush() # flush output. neccessary if stdout ist a pipe (eg | tee). files seem to be ok



def callback_golgi_tendon_organs_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # map from T in the 3D mesh to golgi tendon organs
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_golgi_tendon_organs
  n_output_values = len(output_values[0]) # = n_golgi_tendon_organs (per output slot if there are multiple)
  
  # loop over golgi tendon organs
  for i in range(n_input_values):
    stress = input_values[i]
    
    output_values[0][i] = abs(stress) * 10
    
  # artifical muscle spindle input for debugging
  if False:
    for i in range(n_output_values):
      T = 100   # [ms] cycle duration of stimulus 
      
      if 2000 + i*200 < current_time  <= 2100 + i*200:
        output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 20
  
  print("traction at Golgi tendon organs: {}, output: {}".format(input_values, output_values))

def callback_golgi_tendon_organs_to_interneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping Golgi tendon organs -> interneurons
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_golgi_tendon_organs
  n_output_values = len(output_values[0]) # = n_interneurons (per output slot if there are multiple)
  
  # sum up all input signals and set all outputs to this sum
  
  # collect sum of all Golgi tendon organs
  total_signal = 0
  
  # loop over Golgi tendon organs
  for golgi_tendon_organ_index in range(n_input_values):
    
    # if input is active
    if input_values[golgi_tendon_organ_index] > 20:
      total_signal += input_values[golgi_tendon_organ_index]
    
  # set same value to all connected interneurons
  for interneuron_index in range(n_output_values):
    output_values[0][interneuron_index] = total_signal * 0.01
    
  # initialize buffer the first time
  if 0 not in buffer:
    for golgi_tendon_organ_index in range(n_input_values):
      buffer[golgi_tendon_organ_index] = None
  
  # loop over golgi tendon organ inputs
  for golgi_tendon_organ_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[golgi_tendon_organ_index] > 20:
      buffer[golgi_tendon_organ_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[golgi_tendon_organ_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = 0                   # [ms] delay of the signal
      gaussian_std_dev = 10         # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[golgi_tendon_organ_index]) * 5
      # hodgkin-huxley fires from i_Stim(t) > 4 
        
      output_values[0][golgi_tendon_organ_index] = delayed_signal
  
  print("golgi_tendon_organs_to_interneurons input: {}, output: {}".format(input_values, output_values))

def callback_interneurons_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping interneurons -> input for motor neurons, i.e. signal delay from interneurons to motoneurons
  # the actual N->M mapping from interneurons to motoneurons is done by callback_motoneurons_input
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_interneurons
  n_output_values = len(output_values[0]) # = n_interneurons (per output slot if there are multiple)
  
  # initialize buffer the first time
  if 0 not in buffer:
    for interneuron_index in range(n_input_values):
      buffer[interneuron_index] = None
  
  # loop over interneurons
  for interneuron_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[interneuron_index] > 0:
      buffer[interneuron_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[interneuron_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = golgi_tendon_organ_delay          # [ms] delay of the signal
      gaussian_std_dev = 10                       # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[interneuron_index])   # motor neuron input should be around 1
        
      # loop over output values and set all to the computed signal, cut off at 1e-5
      if delayed_signal > 1e-5:
        #print("interneuron t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[interneuron_index], delayed_signal))
        output_values[0][interneuron_index] = delayed_signal
      else:
        output_values[0][interneuron_index] = None     # signal is below 1e-5, do not set any values
  print("interneurons_to_motoneurons input: {}, output: {}".format(input_values, output_values))



def callback_motoneurons_input(input_values, output_values, current_time, slot_nos, buffer):
    """
    Callback function that transform a number of input_values to a number of output_values.
    This function gets called by a MapDofs object.
    :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
    :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                          The output values buffer, potentially for multiple slots.
                          Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                          the number of required output values. The function should set some of the entries to a computed value.
                          The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
    :param current_time:  Current simulation time.
    :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
    :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                          this function every time. Using this buffer, it is possible to implement a time delay of signals.
    """

    # map from delayed muscle spindle model outputs and delayed interneuron outputs to motoneuron inputs

    # get number of input and output values
    n_input_values = len(input_values)      # = n_muscle_spindles + n_interneurons
    n_output_values = len(output_values[0]) # = n_motoneurons (per output slot if there are multiple)

    # compute mean of all input values
    total_signal = np.mean(input_values)

    # set values to all connected motoneurons
    for motoneuron_index in range(n_output_values):
      # add cortical input as given in the input file, in [nA], motor neuron fires with ~14Hz if drive(t) = 5e-3
      # [timestep, spindle_no]
      timestep_no = int(current_time/1e-2)
      cortical_input_value = cortical_input[timestep_no % cortical_input.shape[0], motoneuron_index % cortical_input.shape[1]]
      
      output_values[0][motoneuron_index] = total_signal + cortical_input_value
      
    print("motoneurons input from spindles and interneurons: {}, resulting drive: {}".format(input_values, output_values))





# TODO, also see get_n_fibers_in_motor_unit
n_motor_units = n_motoneurons


def get_from_obj(data, path):
    for elem in path:
        if type(elem) == str:
            data = data[elem]
        elif type(elem) == int:
            data = data[elem]
        elif type(elem) == tuple:
            # search for key == value with (key, value) = elem
            key, value = elem
            data = next(filter(lambda e: e[key] == value, data))
        else:
            raise KeyError(f"Unknown type of '{elem}': '{type(elem)}'. Path: '{'.'.join(path)}'")
    return data

def muscle1_postprocess(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])
    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    basis_order = get_from_obj(data, [0, 'basisOrder'])
    basis_function = get_from_obj(data, [0, 'basisFunction'])
    assert(basis_function == 'Lagrange')
    assert(basis_order == 2)
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value = 0
    for j in range(ny):
        for i in range(nx):
            z_value += z_data[(nz-1)*nx*ny + j*nx + i]
    z_value /= ny*nx

    global muscle1_tendon_z
    muscle1_tendon_z = z_value
    print("Muscle2: t: {:6.2f}, avg. change of muscle length: {:+2.2f}".format(t, muscle1_tendon_z - muscle1_extent[2]))

