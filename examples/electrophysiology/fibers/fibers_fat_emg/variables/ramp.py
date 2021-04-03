
# scenario name for log file
scenario_name = "ramp"

# material parameters
# --------------------
Conductivity = 3.828      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor units from paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"

import random
random.seed(0)  # ensure that random numbers are the same on every rank
# radius: [μm], stimulation frequency [Hz], jitter [-], activation_start_time [s]
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

end_time = 4000.0                      # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 2.5e-3                        # [ms] timestep width of ODEs (2e-3)
dt_1D = 6.25e-4                        # [ms] timestep width of diffusion (4e-3)
dt_splitting = 2.5e-3                 # [ms] overall timestep width of strang splitting (4e-3)
dt_3D = 5e-1                        # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep_fibers = 50       # [ms] timestep for fiber output, 0.5
output_timestep_3D_emg = 50            # [ms] timestep for output big files of 3D EMG, 100
output_timestep_surface = 10              # [ms] timestep for output surface EMG, 0.5
output_timestep_electrodes = 1    # [ms] timestep for python callback, which is electrode measurement output, has to be >= dt_3D

# input files
fiber_file = "../../../input/left_biceps_brachii_7x7fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_9x9fibers.bin"
#fiber_file = "../../../input/left_biceps_brachii_37x37fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "../../../input/MU_fibre_distribution_10MUs.txt"

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 40      # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40

# HD-EMG electrode parameters
fiber_file_for_hdemg_surface = fat_mesh_file    # use the fat mesh for placing electrodes, this option is the file of the 2D mesh on which electrode positions are set
hdemg_electrode_faces = ["1+"]                  # which faces of this 2D mesh should be considered for placing the HD-EMG electrodes (list of faces, a face is one of "0-" (left), "0+" (right), "1-" (front), "1+" (back))

# xy-direction = across muscle, z-direction = along muscle
hdemg_electrode_offset_xy = 2.0           # [cm] offset from boundary of 2D mesh where the electrode array begins
hdemg_inter_electrode_distance_z = 0.4    # [cm] distance between electrodes ("IED") in z direction (direction along muscle)
hdemg_inter_electrode_distance_xy = 0.4   # [cm] distance between electrodes ("IED") in transverse direction
hdemg_n_electrodes_z = 32           # number of electrodes in z direction (direction along muscle)
hdemg_n_electrodes_xy = 12          # number of electrode across muscle

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = True

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no % len(motor_units)]["radius"]*1e-4
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
