
# scenario name for log file
scenario_name = "heidlauf"

# material parameters
# --------------------
# quantities in mechanics unit system
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)

# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2

c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter

# for debugging, b = 0 leads to normal Mooney-Rivlin
b = 0

material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress

# load
constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
bottom_traction = [0.0,0.0,0.0]        # [N]

# Monodomain parameters
# --------------------
# quantities in CellML unit system
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
sigma_xf = 0                # [mS/cm] conductivity in cross-fiber direction (xf)
sigma_e_f = 6.7             # [mS/cm] conductivity in extracellular space, fiber direction (f)
sigma_e_xf = 3.35           # [mS/cm] conductivity in extracellular space, cross-fiber direction (xf) / transverse

Conductivity = 3.828      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor units from paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
import numpy as np
import random
random.seed(0)  # ensure that random numbers are the same on every rank

n_motor_units = 10   # number of motor units

motor_units = []
for mu_no in range(n_motor_units):

  # capacitance of the membrane
  if mu_no <= 7:
    cm = 0.58    # slow twitch (type I)
  else:
    cm = 1.0     # fast twitch (type II)

  # fiber radius between 40 and 55 [μm]
  min_value = 40
  max_value = 55
  
  # ansatz value(i) = c1 + c2*exp(i),
  # value(0) = min = c1 + c2  =>  c1 = min - c2
  # value(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_value - min_value) / (np.exp(n_motor_units-1) - 1)
  c1 = min_value - c2
  radius = c1 + c2*np.exp(mu_no)

  # stimulation frequency [Hz] between 24 and 7
  min_value = 7 
  max_value = 24
  
  c2 = (max_value - min_value) / (np.exp(n_motor_units-1) - 1)
  c1 = min_value - c2
  stimulation_frequency = c1 + c2*np.exp(n_motor_units-1-mu_no)

  # exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last
  motor_units.append(
  {
    "radius":                radius,                 # [μm] parameter for motor unit: radius of the fiber, used to compute Am
    "cm":                    cm,                     # [uF/cm^2] parameter Cm
    "jitter": [0.1*random.uniform(-1,1) for i in range(100)]     # [-] random jitter values that will be added to the intervals to simulate jitter
  })  

# timing parameters
# -----------------
end_time = 4000.0                      # [ms] end time of the simulation
stimulation_frequency = 1000*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 1e-4                        # [ms] timestep width of ODEs (1e-3)
dt_1D = 1e-4                        # [ms] timestep width of diffusion (1e-3)
dt_splitting = 1e-4                 # [ms] overall timestep width of strang splitting (1e-3)
dt_3D = 1                           # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep_fibers = 4e0       # [ms] timestep for fiber output, 0.5
output_timestep_3D = 4e0              # [ms] timestep for output of fibers and mechanics, should be a multiple of dt_3D


# input files
fiber_file = "cuboid.bin"
firing_times_file = "../../../electrophysiology/input/MU_firing_times_heidlauf_10MU.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "../../../electrophysiology/input/MU_fibre_distribution_10MUs.txt"
cellml_file             = "../../../electrophysiology/input/new_slow_TK_2014_12_08.c"


# fiber mesh
n_fibers_x = 31
n_fibers_y = 13
n_points_whole_fiber = 145  # 600 would be better      # length of muscle is 6 cm

# Heidlauf Diss 5.4.2
# physical size: 6 x 2.9 x (1.2(muscle)+0.2(fat)), number of fibers: 30 x 13, discretization: 144 elements
size_x = 2.9
size_y = 1.2
size_z = 6

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 4      # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40
distribute_nodes_equally = True

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = False
fast_monodomain_solver_optimizations = True # enable the optimizations in the fast multidomain solver
use_analytic_jacobian = True        # If the analytic jacobian should be used for the mechanics problem.
use_vc = True                       # If the vc optimization type should be used for CellmlAdapter

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  #r = motor_units[mu_no % len(motor_units)]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  #return 2./r
  return Am

def get_cm(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["cm"]
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  return stimulation_frequency

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return 0

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return 0
