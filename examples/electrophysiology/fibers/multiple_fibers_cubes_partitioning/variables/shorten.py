
# scenario name for log file
scenario_name = "shorten"

# material parameters
# --------------------

# Fixed units in cellMl models:
# These define the unit system.
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 uA = 1e-6 A
# 1 uF = 1e-6 F
# 
# derived units:
#   (F=s^4*A^2*m^-2*kg^-1) => 1 ms^4*uA^2*cm^-2*x*kg^-1 = (1e-3)^4 s^4 * (1e-6)^2 A^2 * (1e-2)^-2 m^-2 * (x)^-1 kg^-1 = 1e-12 * 1e-12 * 1e4 F = 1e-20 * x^-1 F := 1e-6 F => x = 1e-14
# 1e-14 kg = 10e-15 kg = 10e-12 g = 10 pg

# (N=kg*m*s^-2) => 1 10pg*cm*ms^2 = 1e-14 kg * 1e-2 m * (1e-3)^-2 s^-2 = 1e-14 * 1e-2 * 1e6 N = 1e-10 N = 10 nN
# (S=kg^-1*m^-2*s^3*A^2, Siemens not Sievert!) => (1e-14*kg)^-1*cm^-2*ms^3*uA^2 = (1e-14)^-1 kg^-1 * (1e-2)^-2 m^-2 * (1e-3)^3 s^3 * (1e-6)^2 A^2 = 1e14 * 1e4 * 1e-9 * 1e-12 S = 1e-3 S = 1 mS
# (V=kg*m^2*s^-3*A^-1) => 1 10pg*cm^2*ms^-3*uA^-1 = (1e-14) kg * (1e-2)^2 m^2 * (1e-3)^-3 s^-3 * (1e-6)^-1 A^-1 = 1e-14 * 1e-4 * 1e6 * 1e6 V = 1e-6 V = 1mV
# Hz=s^-1 => 1 ms^-1 = (1e-3)^-1 s^-1 = 1e3 Hz

# Hodgkin-Huxley
# t: ms
# STATES[0], Vm: mV
# CONSTANTS[1], Cm: uF*cm^-2
# CONSTANTS[2], I_Stim: uA*cm^-2
# -> all units are consistent

# Shorten
# t: ms
# CONSTANTS[0], Cm: uF*cm^-2
# STATES[0], Vm: mV
# ALGEBRAIC[32], I_Stim: uA*cm^-2
# -> all units are consistent

sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
sigma_xf = 0                # [mS/cm] conductivity in cross-fiber direction (xf)
sigma_e_f = 6.7             # [mS/cm] conductivity in extracellular space, fiber direction (f)
sigma_e_xf = 3.35           # [mS/cm] conductivity in extracellular space, cross-fiber direction (xf) / transverse

Conductivity = sigma_f              # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

import random
motor_units = [
  {"radius": 40.00, "activation_start_time": 0.0, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # low number of fibers
  {"radius": 42.35, "activation_start_time": 0.2, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 45.00, "activation_start_time": 0.4, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 48.00, "activation_start_time": 0.6, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 51.42, "activation_start_time": 0.8, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 55.38, "activation_start_time": 1.0, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 60.00, "activation_start_time": 1.2, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 65.45, "activation_start_time": 1.4, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 72.00, "activation_start_time": 1.6, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 80.00, "activation_start_time": 1.8, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # high number of fibers
]

# timing parameters
# -----------------
end_time = 4000.0                      # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
output_timestep = 3               # [ms] timestep for output files, 0.5
dt_0D = 1e-3                        # [ms] timestep width of ODEs
dt_1D = 1e-3                        # [ms] timestep width of diffusion
dt_splitting = 1e-3                 # [ms] overall timestep width of strang splitting

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
#fiber_file = "../../../input/left_biceps_brachii_1x1fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_9x9fibers.bin"
firing_times_file="../../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file="../../../input/MU_fibre_distribution_10MUs.txt"
cellml_file = "../../../input/new_slow_TK_2014_12_08.c"

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r
  #return Am

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
  return 0
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
