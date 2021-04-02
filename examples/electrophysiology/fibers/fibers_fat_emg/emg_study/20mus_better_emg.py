import os

# scenario name for log file
scenario_name = "20mus"

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
import numpy as np

n_motor_units = 20   # number of motor units

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

  # stimulation frequency [Hz] between 24 and 7
  min_value = 7
  max_value = 24
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  stimulation_frequency = c1 + c2*1.02**(n_motor_units-1-mu_no)

  # exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last
  motor_units.append(
  {
    "radius":                radius,                 # [μm] parameter for motor unit: radius of the fiber, used to compute Am
    "cm":                    cm,                     # [uF/cm^2] parameter Cm
    "activation_start_time": 1*mu_no,              # [s] when to start activating this motor unit, here it is a ramp
    "stimulation_frequency": stimulation_frequency,  # [Hz] stimulation frequency for activation
    "jitter": [0.1*random.uniform(-1,1) for i in range(100)]     # [-] random jitter values that will be added to the intervals to simulate jitter
  })  

# timing parameters
# -----------------
end_time = 40000.0                      # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 2.5e-3                        # [ms] timestep width of ODEs (2e-3)
dt_1D = 6.25e-4                        # [ms] timestep width of diffusion (4e-3)
dt_splitting = 2.5e-3                 # [ms] overall timestep width of strang splitting (4e-3)
dt_3D = 5e-1                        # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep_fibers = 2e5       # [ms] timestep for fiber output, 0.5
output_timestep_3D_emg = 2e5            # [ms] timestep for output big files of 3D EMG, 100
output_timestep_surface = 1e5              # [ms] timestep for output surface EMG, 0.5
output_timestep_electrodes = 2e8    # [ms] timestep for python callback, which is electrode measurement output, has to be >= dt_3D

# input files
input_directory = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")
fiber_file              = input_directory+"/left_biceps_brachii_37x37fibers.bin"
fat_mesh_file           = input_directory+"/left_biceps_brachii_37x37fibers_thin_fat.bin"
firing_times_file       = input_directory+"/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "MU_fibre_distribution_37x37_20a.txt"

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
disable_firing_output = False
fast_monodomain_solver_optimizations = True # enable the optimizations in the fast multidomain solver
use_vc = True                       # If the vc optimization type should be used for CellmlAdapter

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no % len(motor_units)]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r

def get_cm(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["cm"]
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  #return 1
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3

# output motor_unit config in a readable format
if True:
  import sys
  # only on rank 0
  if (int)(sys.argv[-2]) == 0:
    for mu_no,item in enumerate(motor_units[0::4]):
      print("MU {}".format(mu_no*4))
      for (key,value) in item.items():
        if key != "jitter":
          print("  {}: {}".format(key,value))
