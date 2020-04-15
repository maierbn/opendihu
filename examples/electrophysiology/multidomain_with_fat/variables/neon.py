# scenario name for log file
scenario_name = "neon"

# material parameters
# --------------------
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
sigma_xf = 0                # [mS/cm] conductivity in cross-fiber direction (xf)
sigma_e_f = 6.7             # [mS/cm] conductivity in extracellular space, fiber direction (f)
sigma_e_xf = 3.35           # [mS/cm] conductivity in extracellular space, cross-fiber direction (xf) / transverse

Am = 500.0                  # [cm^-1] surface area to volume ratio
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)

# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor units from paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
import random
random.seed(0)  # ensure that random numbers are the same on every rank
#   fiber_no: center MU around this fiber
#   standard_deviation [-]: relative to muscle diameter, 
#   maximum [-]: create f_r as gaussian with standard_deviation and maximum around the fiber given in fiber_no
#   radius: [μm], activation_start_time: [s], stimulation frequency [Hz], jitter [-]
motor_units = [
  {"fiber_no": 10, "standard_deviation": 0.2, "maximum": 0.2, "radius": 40.00, "activation_start_time": 0.0, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # low number of fibers
  {"fiber_no": 20, "standard_deviation": 0.2, "maximum": 0.2, "radius": 42.35, "activation_start_time": 0.2, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 30, "standard_deviation": 0.2, "maximum": 0.2, "radius": 45.00, "activation_start_time": 0.4, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 40, "standard_deviation": 0.2, "maximum": 0.2, "radius": 48.00, "activation_start_time": 0.6, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 55, "standard_deviation": 0.2, "maximum": 0.2, "radius": 51.42, "activation_start_time": 0.8, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 60, "standard_deviation": 0.2, "maximum": 0.2, "radius": 55.38, "activation_start_time": 1.0, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 70, "standard_deviation": 0.2, "maximum": 0.2, "radius": 60.00, "activation_start_time": 1.2, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 80, "standard_deviation": 0.2, "maximum": 0.2, "radius": 65.45, "activation_start_time": 1.4, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 50, "standard_deviation": 0.2, "maximum": 0.2, "radius": 72.00, "activation_start_time": 1.6, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"fiber_no": 25, "standard_deviation": 0.2, "maximum": 0.2, "radius": 80.00, "activation_start_time": 1.8, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # high number of fibers
]
# solvers
# -------
multidomain_solver_type = "gmres"          # solver for the multidomain problem
multidomain_preconditioner_type = "ilu"   # preconditioner

# set initial guess to zero for direct solver
initial_guess_nonzero = "lu" not in multidomain_solver_type 


# timing parameters
# -----------------
end_time = 4000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 3e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 3e-3               # [ms] timestep width of the multidomain solver
dt_splitting = 3e-3                 # [ms] overall timestep width of strang splitting (3e-3)
output_timestep_multidomain = 2e-1  # [ms] timestep for multidomain output
output_timestep_multidomain = 1     # [ms] timestep for multidomain output
#end_time = 1e-2

# input files
cellml_file = "../../input/hodgkin_huxley_1952.c"
fiber_file = "../../input/left_biceps_brachii_9x9fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../input/MU_firing_times_immediately.txt"


# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 20
sampling_stride_fat = 1

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = False

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r
  #return Am

def get_cm(mu_no):
  return Cm
  
def get_conductivity(mu_no):
  return Conductivity

def get_specific_states_call_frequency(mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(mu_no):
  #return 0
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
