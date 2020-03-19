
# scenario name for log file
scenario_name = "debug"

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
# motor unit parameters
# stimulation frequency in [Hz], activation_start_time in [s]
#   fiber_no: center MU around this fiber, create f_r as gaussion from standard_deviation and maximum

motor_units = [
  {"fiber_no": 10, "standard_deviation": 20.0, "maximum": 0.5, "activation_start_time": 0.0, "stimulation_frequency": 10.0,},
  {"fiber_no": 30, "standard_deviation": 20.0, "maximum": 0.4, "activation_start_time": 0.0, "stimulation_frequency": 10.0,},
  {"fiber_no": 40, "standard_deviation": 30.0, "maximum": 0.6, "activation_start_time": 0.0, "stimulation_frequency": 10.0,},
]
motor_units = motor_units[0:2]

end_time = 4000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 3e-3                        # [ms] timestep width of ODEs (2e-3)
dt_splitting = 3e-3                 # [ms] overall timestep width of strang splitting (4e-3)
output_timestep = 2e-1              # [ms] timestep for output big files of 3D EMG, 100
output_timestep = 1e-1              # [ms] timestep for output big files of 3D EMG, 100

# input files
cellml_file = "../../input/hodgkin_huxley_1952.c"
fiber_file = "../../input/left_biceps_brachii_7x7fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../input/MU_firing_times_immediately.txt"


# stride for sampling the 3D elements from the fiber data
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 20
#sampling_stride_z = 50   # faster, but stimulus does not propagate

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = False

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(mu_no):
  return Am

def get_cm(mu_no):
  return Cm
  
def get_conductivity(mu_no):
  return Conductivity

def get_specific_states_call_frequency(mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(mu_no):
  return 0
  #return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
