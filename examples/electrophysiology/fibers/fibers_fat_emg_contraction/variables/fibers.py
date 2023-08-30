
# scenario name for log file
scenario_name = "fibers"

# Fixed units in cellMl models:
# These define the unit system.
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 uA = 1e-6 A
# 1 uF = 1e-6 F
# 
# Fixed units in mechanics system
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 N
# 1 N/cm^2 = (kg*m*s^-2) / (1e-2 m)^2 = 1e4 kg*m^-1*s^-2 = 10 kPa
# (kg = N*s^2*m^-1) => N*ms^2*cm^-1 = N*(1e-3 s)^2 * (1e-2 m)^-1 = 1e-4 N*s^2*m^-1 = 1e-4 kg
# (kg/m^3) => 1 * 1e-4 kg * (1e-2 m)^-3 = 1e2 kg/m^3
# (m/s^2) => 1 cm/ms^2 = 1e-2 m * (1e-3 s)^-2 = 1e4 m*s^-2

# material parameters
# --------------------
# quantities in mechanics unit system

# parameters for main simulation
# load
main_constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
main_bottom_traction = [0,0,-1e-3]        # [N]  (-30 works)

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
#b = 0

material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress
#Conductivity = 3.828        # [mS/cm] sigma, conductivity 
Conductivity = 8.93         # [mS/cm] sigma, conductivity 

# timing and activation parameters
# -----------------
# motor unit parameters similar to paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
# however, values from paper fail for mu >= 6, then stimulus gets reflected at the ends of the muscle, therefore fiber radius is set to <= 55

import random
random.seed(0)  # ensure that random numbers are the same on every rank
import numpy as np

n_fibers_in_fiber_file = 81
n_motor_units = 10   # number of motor units

motor_units = []
for mu_no in range(n_motor_units):

  # capacitance of the membrane
  if mu_no <= 0.7*n_motor_units:
    cm = 0.58    # slow twitch (type I)
  else:
    cm = 1.0     # fast twitch (type II)

  # fiber radius between 40 and 55 [μm]
  min_value = 40
  max_value = 55

  # ansatz value(i) = c1 + c2*exp(i),
  # value(0) = min = c1 + c2  =>  c1 = min - c2
  # value(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  radius = c1 + c2*1.02**(mu_no)

  # standard_deviation
  min_value = 0.1
  max_value = 0.6
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  standard_deviation = c1 + c2*1.02**mu_no
  maximum = 2.0/n_motor_units*standard_deviation

  # exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last
  motor_units.append(
  {
    "fiber_no":              random.randint(0,n_fibers_in_fiber_file),  # [-] fiber from input files that is the center of the motor unit domain
    "maximum":               maximum,                # [-] maximum value of f_r, create f_r as gaussian with standard_deviation and maximum around the fiber 
    "standard_deviation":    standard_deviation,     # [-] standard deviation of f_r
    "radius":                radius,                 # [μm] parameter for motor unit: radius of the fiber, used to compute Am
    "cm":                    cm,                     # [uF/cm^2] parameter Cm
  })

#motor_units = motor_units[0:1]  # for debugging, only 1 motor unit

# solvers
# -------
# 3D potential flow
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "none" # preconditioner

# 3D EMG, bidomain solver
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False   # If the initial guess for the emg linear system should be set to the previous solution
emg_solver_maxit = 1e4              # maximum number of iterations
emg_solver_abstol = 1e-5            # absolute tolerance of the residual of the linear solver
emg_solver_reltol = 1e-5            # relative tolerance of the residual of the linear solver

# elasticity
elasticity_solver_type = "lu"
elasticity_preconditioner_type = "none"
snes_max_iterations = 34                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 1       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-4      # absolute tolerance of the nonlinear solver
relative_tolerance = 1e-10           # relative tolerance of the residual of the linear solver
absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver

theta = 1.0                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True   # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations

# timing parameters
# -----------------
end_time = 5_000.0                  # [ms] end time of the simulation
dt_0D = 2.5e-5                        # [ms] timestep width of ODEs (1e-3), for shorten use 2.5e-5
dt_1D = 2.5e-5                      # [ms] timestep width of the 1D electric conduction problem, for shorten use 2.5e-5
dt_splitting = 2.5e-5                # [ms] timestep width of strang splitting between 0D and 1D for the fibers, for shorten use 2.5e-5
dt_3D = 1e-3                        # [ms] timestep width of the bidomain solver
dt_elasticity = 1e-1                # [ms] time step width of elasticity solver
#dt_elasticity = 1e-2                # [ms] time step width of elasticity solver

dt_neurons = 1e-2                   # [ms] same timestep width for all neuron solvers
dt_golgi_tendon_organs = dt_neurons # [ms] timestep width of cellml solver of golgi tendon organs
dt_muscle_spindles     = 1e-3       # [ms] timestep width of cellml solver of muscle spindles
dt_interneuron         = dt_neurons # [ms] timestep width of the cellml solver for interneurons
dt_motoneuron          = dt_neurons # [ms] timestep width of the cellml solver for motoneurons

dt_neuron_transfer     = dt_elasticity  # [ms] interval when to call callback functions and transfer values between CellML models, increase this to speed up the simulation
#dt_neuron_transfer     = dt_neurons  # [ms] interval when to call callback functions and transfer values between CellML models, increase this to speed up the simulation

output_timestep_elasticity = 2      # [ms] timestep for elasticity output files
output_timestep_fibers = 2         # [ms] timestep for fiber output files
output_timestep_3D_emg = 2         # [ms] timestep for output of 3D emg

#output_timestep_multidomain = dt_elasticity
#output_timestep_elasticity = dt_elasticity

# input files
# -----------
import os
import opendihu
input_directory   = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")

# if the hodgkin_huxley-razumova model is used
if "hh" in opendihu.program_name:
  cellml_file       = input_directory+"/hodgkin_huxley-razumova.cellml"
  dt_0D = 1e-3                   # [ms] timestep width of ODEs (1e-3), for shorten use 2.5e-5
  dt_1D = 1e-3                   # [ms] timestep width of the 1D electric conduction problem, for shorten use 2.5e-5
  dt_splitting = 1e-3            # [ms] timestep width of strang splitting between 0D and 1D for the fibers, for shorten use 2.5e-5
  scenario_name = "spindles_fibers_hh"

else:
  cellml_file       = input_directory+"/new_slow_TK_2014_12_08.c"
  scenario_name = "fibers"

fiber_file        = input_directory+"/left_biceps_brachii_9x9fibers_b.bin"  # this is a variant of 9x9fibers with a slightly different mesh that somehow works better
fiber_file        = input_directory+"/left_biceps_brachii_9x9fibers.bin"
fat_mesh_file     = fiber_file + "_fat.bin"

fiber_file = "cuboid.bin"
fat_mesh_file = "cuboid_fat2.bin"


# fiber mesh
n_fibers_x = 31
n_fibers_y = 13
n_points_whole_fiber = 145  # 600 would be better      # length of muscle is 6 cm

# Heidlauf Diss 5.4.2
# physical size: 6 x 2.9 x (1.2(muscle)+0.2(fat)), number of fibers: 30 x 13, discretization: 144 elements
size_x = 2.9
size_y = 1.2
size_z = 6


firing_times_file = input_directory+"/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
firing_times_file = input_directory+"/MU_firing_times_once.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
firing_times_file = input_directory+"/MU_firing_times_heidlauf_10MU.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = input_directory+"/MU_fibre_distribution_10MUs.txt"
cortical_input_file = input_directory+"/cortical_input_realistic.txt"

# comment to cortical_input_realistic.txt:
# - mean cortical drive: 8nA
# - common noise: 15-35Hz filtered white noise, Amplitude 16% CoV (relativ with respect to the mean value).
# - independent noise: bandpass filtered white noise (100Hz low pass), standart deviation 25% (relativ to the common noise).  
# Jede Zeile repräsentiert einen Zeitschritt (0.01ms, d.h. Abtastrate 100000Hz, insgesamt sollten es dann 10s Signal sein) und jede Spalte repräsentiert eine motorische Einheit (N_MU=10).


# stride for meshes
# -----------------

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
# If you change this, delete the compartment_relative_factors.* files, they have to be generated again.
sampling_stride_x = 1 
sampling_stride_y = 1 
sampling_stride_z = 10
sampling_stride_fat = 1 

# how much of the 3D mesh is used for elasticity
sampling_factor_elasticity_x = 0.7 
sampling_factor_elasticity_y = 0.7 
sampling_factor_elasticity_z = 0.7 
sampling_factor_elasticity_fat_y = 0.5 

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
states_output = False                    # if also the subcellular states should be output, this produces large files, set output_timestep_0D_states
enable_surface_emg = True               # if the EMG values on a 2D surface should be written to files
show_linear_solver_output = False       # if every solve of multidomain diffusion should be printed
disable_firing_output = False            # if information about firing of MUs should be printed
optimization_type = "vc"                # the optimization_type used in the cellml adapter, "vc" uses explicit vectorization
approximate_exponential_function = True # if the exponential function should be approximated by a Taylor series with only 11 FLOPS

# fiber callbacks
# ----------------------
# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
# the diffusion prefactor is computed by get_diffusion_prefactor at helper.py:764  as conductivity/(am*cm)

def get_am(fiber_no, mu_no):
  if mu_no >= len(motor_units) or mu_no >= n_motor_units:
    print("Warning, mu_no={}, n_motor_units={}={}".format(mu_no,n_motor_units,len(motor_units)))
    
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no % len(motor_units)]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r

def get_cm(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["cm"]
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity
  

def get_specific_states_call_frequency(fiber_no, mu_no):
  return 1
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return 0
  #return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  #return 1000  # start directly
  return 0  # start directly
  #return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
