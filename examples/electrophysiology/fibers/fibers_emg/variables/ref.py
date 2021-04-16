
# scenario name for log file
scenario_name = "ramp_emg"

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

PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]

# timing parameters
# -----------------
end_time = 4000.0                      # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 3e-3                        # [ms] timestep width of ODEs (2e-3)
dt_1D = 1e-3                        # [ms] timestep width of diffusion (4e-3)
dt_splitting = 3e-3                 # [ms] overall timestep width of strang splitting (4e-3)
dt_3D = 4e-1                        # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep = 4e-1              # [ms] timestep for output files, 0.5
output_timestep_fibers = 4e-1       # [ms] timestep for output files, 0.5
output_timestep_3D_emg = 25            # [ms] timestep for output big files, 100

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
#fiber_file="../../../input/7x7fibers.bin"
#fiber_file="../../../input/37x37fibers.bin"
#fiber_file="../../../input/left_biceps_brachii_37x37fibers.bin"
#fiber_file = "../../../input/left_biceps_brachii_7x7fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_7x7fibers.bin"
firing_times_file="../../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file="../../../input/MU_fibre_distribution_10MUs.txt"

fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../../input/MU_firing_times_immediately.txt"

cellml_file = "../../../input/hodgkin_huxley_1952.c"

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  return Am

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  return stimulation_frequency

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return 0

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return 0
