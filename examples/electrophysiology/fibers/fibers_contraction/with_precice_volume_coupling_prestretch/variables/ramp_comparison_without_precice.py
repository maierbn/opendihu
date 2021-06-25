
# scenario name for log file and output directory under "out"
scenario_name = "ramp"

# material parameters
# --------------------
# parameters for precontraction
# -----------------------------
# load
#precontraction_constant_body_force = (0,0,50*9.81e-4)   # [cm/ms^2], gravity constant for the body force
precontraction_constant_body_force = (0,0,10*9.81e-4)   # [cm/ms^2], gravity constant for the body force
precontraction_bottom_traction = [0,0,0]        # [N]
constant_gamma = 0.25    # if it diverges, reduce this value, the active stress will be pmax*constant_gamma

# parameters for prestretch
# -----------------------------
# load
prestretch_constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
prestretch_bottom_traction = [0,0,-10]        # [N]  (-30 also works)

# parameters for main simulation
# load
main_constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
main_bottom_traction = [0,0,-10]        # [N]  (-30 works)

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
material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress

# for debugging, b = 0 leads to normal Mooney-Rivlin
#b = 0

material_parameters = [c1, c2, b, d]   # material parameters

# quantities in CellML unit system
#Conductivity = 3.828        # [mS/cm] sigma, conductivity 
Conductivity = 8.93         # [mS/cm] sigma, conductivity 
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)

Conductivity = sigma_f      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force

# timing and activation parameters
# -----------------
import random
random.seed(0)  # ensure that random numbers are the same on every rank
motor_units = [
  {"radius": 40.00, "activation_start_time": 0.00, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #1    # low number of fibers
  {"radius": 42.35, "activation_start_time": 0.02, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #2
  {"radius": 45.00, "activation_start_time": 0.04, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #3
  {"radius": 48.00, "activation_start_time": 0.06, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #4
  {"radius": 51.42, "activation_start_time": 0.08, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #5
  {"radius": 55.38, "activation_start_time": 0.10, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #6 
  {"radius": 60.00, "activation_start_time": 0.12, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #7
  {"radius": 65.45, "activation_start_time": 0.14, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #8
  {"radius": 72.00, "activation_start_time": 0.16, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #9
  {"radius": 80.00, "activation_start_time": 0.18, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #10   # high number of fibers
]
# note: negative start time is the same as zero, it is just there for debugging. Delete the minus signs to get a ramp

end_time = 10                      # [ms] end time of the simulation
dt_0D = 2.5e-5                      # [ms] timestep width of ODEs (2e-3)
dt_1D = 2.5e-5                      # [ms] timestep width of diffusion (4e-3)
dt_splitting = 2.5e-5               # [ms] overall timestep width of strang splitting (4e-3)
dt_3D = 1                           # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG, this has to be the same value as in the precice_config.xml
output_timestep = 10.0                # [ms] timestep for output files, 5.0
output_timestep_fibers = output_timestep   # [ms] timestep for fiber output
output_timestep_big = 100.0            # [ms] timestep for output big files of 3D EMG, 100
output_timestep_elasticity = 20      # [ms] timestep for elasticity output files

# The values of dt_3D and end_time have to be also defined in "precice-config.xml" with the same value (the value is only significant in the precice-config.xml, the value here is used for output writer time intervals)
# <max-time value="100.0"/>           <!-- end time of the whole simulation -->
# <time-window-size value="1e0"/>   <!-- timestep width dt_3D -->

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements

import opendihu

# parameters for the contraction program
if "contraction" in opendihu.program_name:
  sampling_stride_x = 2
  sampling_stride_y = 2
  sampling_stride_z = 185
  # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40 

else:
  # parameters for the fibers_with_3d program
  
  sampling_stride_x = 2
  sampling_stride_y = 2
  sampling_stride_z = 185
  # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40 

distribute_nodes_equally = True     # (default: False)
# True: set high priority to make subdomains have approximately equal number of fibers but creates tiny remainder elements inside the subdomains
# False: make elements more equally sized, this can lead to a slight imbalance in the number of fibers per subdomain

# Tolerance value in the element coordinate system of the 3D elements, [0,1]^3
# when a fiber point is still considered part of the element.
# Try to increase this such that all mappings have all points.
mapping_tolerance = 0.3
#mapping_tolerance = 1.0

# input files
# -----------
import os
import opendihu
input_directory   = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")

#fiber_file        = input_directory + "/left_biceps_brachii_7x7fibers.bin"
fiber_file        = input_directory + "/left_biceps_brachii_9x9fibers.old.bin"
#fiber_file        = input_directory + "/left_biceps_brachii_13x13fibers.bin"
firing_times_file = input_directory + "/MU_firing_times_always.txt"
#fiber_distribution_file = input_directory + "/MU_fibre_distribution_10MUs.txt"
fiber_distribution_file = input_directory + "/MU_fibre_distribution_9x9_10.txt"
#fiber_distribution_file = input_directory + "/MU_fibre_distribution_10MUs_13x13.txt"

# if the hodgkin_huxley-razumova model is used
if "hh" in opendihu.program_name:
  cellml_file       = input_directory+"/hodgkin_huxley-razumova.cellml"
  dt_0D = 1e-3                   # [ms] timestep width of ODEs (1e-3), for shorten use 2.5e-5
  dt_1D = 1e-3                   # [ms] timestep width of the 1D electric conduction problem, for shorten use 2.5e-5
  dt_splitting = 1e-3            # [ms] timestep width of strang splitting between 0D and 1D for the fibers, for shorten use 2.5e-5
  scenario_name = "ramp"

else:
  cellml_file       = input_directory+"/new_slow_TK_2014_12_08.c"
  scenario_name = "ramp"

# EMG solver parameters
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = True    # If the initial guess for the emg linear system should be set to the previous solution
emg_solver_maxit = 1e4              # maximum number of iterations for the static bidomain solver
emg_solver_reltol = 1e-5            # relative tolerance for solver
emg_solver_abstol = 1e-5            # absolute tolerance for solver

# other options
states_output = False           # if the 0D states should be written (this produces large output files)
paraview_output = True         # produce output files for paraview
adios_output = False           # produce adios output files
exfile_output = False          # produce exfiles output files
python_output = False          # produce python output files


# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3

# callback function for artifical stress values in precontraction computation
def set_gamma_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, custom_argument):
    # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
    # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
    #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
    # time_step_no:        (int)   current time step number
    # current_time:        (float) the current simulation time
    # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
    #                       i.e. [point0_component0, point0_component1, ... point0_componentN, point1_component0, point1_component1, ...]
    #                       After the call, these values will be assigned to the field variable.
    # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
    # additional_argument: The value of the option "additionalArgument", can be any Python object.
    
    # set all values to 1
    for i in range(len(values)):
      values[i] = constant_gamma

# callback function for artifical lambda values in precontraction computation
def set_lambda_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, custom_argument):
    # set all values to 1
    for i in range(len(values)):
      values[i] = 1.0
