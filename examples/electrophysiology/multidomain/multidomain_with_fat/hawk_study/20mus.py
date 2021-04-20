# this is a combined variables file that works for fibers_with_fat_emg and multidomain_with_fat
# You have to specify multidomain or fibers as third argument, e.g.:
#   multidomain_with_fat_emg ../settings_multidomain_with_fat.py 20mus.py multidomain
#   fibers_fat_emg ../settings_fibers_fat_emg.py 20mus.py fibers

import sys,os,time

# import the script generate_fibers_analogous_to_multidomain
sys.path.insert(0, os.path.dirname(__file__))
import generate_fibers_analogous_to_multidomain

# parse command line arguments
is_multidomain = "multidomain" in sys.argv[1]
rank_no = (int)(sys.argv[-2])

if rank_no == 0:
  print("This simulation is with {}.".format("multidomain" if is_multidomain else "fibers"))

# scenario name for log file
scenario_name = "20mus_multidomain" if is_multidomain else "20mus_fibers" 

# material parameters
# --------------------
Conductivity = 3.828      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor unit parameters similar to paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
# however, values from paper fail for mu >= 6, then stimulus gets reflected at the ends of the muscle, therefore fiber radius is set to <= 55

import random
random.seed(0)  # ensure that random numbers are the same on every rank
import numpy as np

n_fibers_x = 31
if not is_multidomain:   # for fibers use finer mesh here
  n_fibers_x = 31
n_fibers_in_fiber_file = n_fibers_x*n_fibers_x
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
  maximum = 10.0/n_motor_units*standard_deviation

  # stimulation frequency [Hz] between 24 and 7
  min_value = 7
  max_value = 24
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  stimulation_frequency = c1 + c2*1.02**(n_motor_units-1-mu_no)

  # exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last
  motor_units.append(
  {
    "fiber_no":              random.randint(0,n_fibers_in_fiber_file),  # [-] fiber from input files that is the center of the motor unit domain
    "maximum":               maximum,                # [-] maximum value of f_r, create f_r as gaussian with standard_deviation and maximum around the fiber g
    "standard_deviation":    standard_deviation,     # [-] standard deviation of f_r
    "radius":                radius,                 # [μm] parameter for motor unit: radius of the fiber, used to compute Am
    "cm":                    cm,                     # [uF/cm^2] parameter Cm
    "activation_start_time": 1e-3*10*mu_no,             # [s] when to start activating this motor unit, here it is a ramp
    "stimulation_frequency": stimulation_frequency,  # [Hz] stimulation frequency for activation
    "jitter": [0.1*random.uniform(-1,1) for i in range(100)]     # [-] random jitter values that will be added to the intervals to simulate jitter
  })

# rescale maximum such that the total sum of the probabilities is equal to 1
maximum_probability_sum = generate_fibers_analogous_to_multidomain.get_maximum_probability_sum(n_fibers_x, motor_units)

if rank_no == 0:
  print("divide maximum f_r of motor units by {}".format(maximum_probability_sum))
for i in range(len(motor_units)):
  motor_units[i]["maximum"] /= maximum_probability_sum


# solvers
# -------
# multidomain
multidomain_solver_type = "gmres"          # solver for the multidomain problem
#multidomain_solver_type = "lu"          # solver for the multidomain problem
#multidomain_preconditioner_type = "bjacobi"   # preconditioner
#multidomain_preconditioner_type = "boomeramg"   # preconditioner
multidomain_preconditioner_type = "euclid"   # euclid ilu preconditioner, has a memory leak
#multidomain_preconditioner_type = "none"   # sor
#multidomain_preconditioner_type = "sor"   # sor

multidomain_alternative_solver_type = "gmres"            # alternative solver, used when normal solver diverges
multidomain_alternative_preconditioner_type = "euclid"    # preconditioner of the alternative solver

# set initial guess to zero for direct solver
initial_guess_nonzero = "lu" not in multidomain_solver_type 

# if using boomeramg, tolerances cannot be as low as 1e-10, otherwise it becomes unstable
multidomain_absolute_tolerance = 1e-15    # absolute residual tolerance for the multidomain solver
multidomain_relative_tolerance = 1e-15    # relative residual tolerance for the multidomain solver
theta = 1.0                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True   # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
multidomain_max_iterations = 1e3                         # maximum number of iterations
multidomain_alternative_solver_max_iterations = 1e4      # maximum number of iterations of the alternative solver

# timing parameters
# -----------------
end_time = 1000.0                   # [ms] end time of the simulation
end_time = 200
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.

if is_multidomain:
  # multidomain
  dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
  dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver
  dt_splitting = 1e-3                 # [ms] overall timestep width of strang splitting (3e-3)
  output_timestep_multidomain = 0.5     # [ms] (10) timestep for multidomain output
  #output_timestep_0D_states = 2e-3    # [ms] timestep for output files of 0D subcellular model states

else:
  # fibers
  dt_0D = 2.5e-3                        # [ms] timestep width of ODEs (2e-3)
  dt_1D = 6.25e-4                        # [ms] timestep width of diffusion (4e-3)
  dt_splitting = 2.5e-3                 # [ms] overall timestep width of strang splitting (4e-3)
  dt_3D = 5e-1                        # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
  output_timestep_fibers = 1e2       # [ms] timestep for fiber output, 0.5
  output_timestep_3D_emg = 1e2            # [ms] timestep for output big files of 3D EMG, 100
  output_timestep_electrodes = 2e8    # [ms] timestep for python callback, which is electrode measurement output, has to be >= dt_3D

output_timestep_surface = 2              # [ms] timestep for output surface EMG, 0.5

# input files
input_directory = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")
cellml_file             = input_directory + "/hodgkin_huxley_1952.c"
#fiber_file              = input_directory + "/left_biceps_brachii_13x13fibers.bin"
fiber_file              = input_directory + "/left_biceps_brachii_{}x{}fibers.bin".format(n_fibers_x, n_fibers_x)
fat_mesh_file           = fiber_file + "_fat.bin"
firing_times_file       = input_directory + "/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency

# generate fiber distribution file
fiber_distribution_file = "MU_fiber_distribution_{}x{}_{}mus.txt".format(n_fibers_x,n_fibers_x,n_motor_units)
if os.path.exists(fiber_distribution_file):
  if rank_no == 0:
    print("Fiber distribution file \"{}\" already exists, do not create again.".format(fiber_distribution_file))
else:
  if rank_no == 0:
    print("Create fiber distribution file \"{}\".".format(fiber_distribution_file))
    generate_fibers_analogous_to_multidomain.generate_fiber_file(fiber_distribution_file, n_fibers_x, motor_units)
  else:
    # wait while the file does not exists
    while not os.path.exists(fiber_distribution_file):
      time.sleep(1)


# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 20      # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40
local_sampling_stride_z = 1 # stride within the local domain, leads to potentially unequal element sizes in z direction if set != 1
sampling_stride_fat = 1

# HD-EMG electrode parameters
fiber_file_for_hdemg_surface = fat_mesh_file    # use the fat mesh for placing electrodes, this option is the file of the 2D mesh on which electrode positions are set
hdemg_electrode_faces = ["1+"]                  # which faces of this 2D mesh should be considered for placing the HD-EMG electrodes (list of faces, a face is one of "0-" (left), "0+" (rig

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
states_output = False    # if also the subcellular states should be output, this produces large files, set output_timestep_0D_states
show_linear_solver_output = False    # if every solve of multidomain diffusion should be printed
disable_firing_output = False   # if information about firing of MUs should be printed
fast_monodomain_solver_optimizations = True # enable the optimizations in the fast multidomain solver
use_vc = True                       # If the vc optimization type should be used for CellmlAdapter
debug_output = False     # enable additional debugging output during partitioning


if is_multidomain:
  # functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
  def get_am(mu_no):
    return Am
    # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
    r = motor_units[mu_no % len(motor_units)]["radius"]*1e-4
    # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
    return 2./r

  def get_cm(mu_no):
    return Cm
    return motor_units[mu_no % len(motor_units)]["cm"]
    
  def get_conductivity(mu_no):
    return Conductivity

  def get_specific_states_call_frequency(mu_no):
    stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
    return stimulation_frequency*1e-3

  def get_specific_states_frequency_jitter(mu_no):
    return motor_units[mu_no % len(motor_units)]["jitter"]

  def get_specific_states_call_enable_begin(mu_no):
    #return 0  # start directly
    return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
else:
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
    print("")
    for mu_no,item in enumerate(motor_units):    
      print("MU {:2d} around fiber {:3d}, fr_max: {:.3f}, stddev: {:.2f}, r: {:.1f} μm, Cm: {:.1f} uF/cm^2, tstart: {:.3f} s, f: {:.3f} Hz".
        format(mu_no, item["fiber_no"], item["maximum"], item["standard_deviation"], item["radius"], item["cm"], 
               item["activation_start_time"], item["stimulation_frequency"]))
    print("")
