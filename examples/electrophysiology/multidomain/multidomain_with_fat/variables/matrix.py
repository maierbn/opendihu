# scenario name for log file
scenario_name = "matrix"

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

n_fibers_in_fiber_file = 169
n_motor_units = 3   # number of motor units

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
    "activation_start_time": 0.1*mu_no,              # [s] when to start activating this motor unit, here it is a ramp
    "stimulation_frequency": stimulation_frequency,  # [Hz] stimulation frequency for activation
    "jitter": [0.1*random.uniform(-1,1) for i in range(100)]     # [-] random jitter values that will be added to the intervals to simulate jitter
  })


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
end_time = 2e-3                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver
dt_splitting = 1e-3                 # [ms] overall timestep width of strang splitting (3e-3)
output_timestep_multidomain = 1e-3     # [ms] timestep for multidomain output
#output_timestep_0D_states = 2e-3    # [ms] timestep for output files of 0D subcellular model states


# input files
cellml_file = "../../../input/hodgkin_huxley_1952.c"
fiber_file = "../../../input/left_biceps_brachii_9x9fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_13x13fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../../input/MU_firing_times_immediately.txt"
firing_times_file = "../../../input/MU_firing_times_always.txt"


# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 3
sampling_stride_y = 3
sampling_stride_z = 20      # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40
sampling_stride_fat = 1

distribute_nodes_equally = True     # (default: False)
# True: set high priority to make subdomains have approximately equal number of fibers but creates tiny remainder elements inside the subdomains
# False: make elements more equally sized, this can lead to a slight imbalance in the number of fibers per subdomain

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

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no % len(motor_units)]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r
  #return Am

def get_cm(mu_no):
  return motor_units[mu_no % len(motor_units)]["cm"]
  #return Cm
  
def get_conductivity(mu_no):
  return Conductivity

def get_specific_states_call_frequency(mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(mu_no):
  #return 0
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(mu_no):
  #return 0  # start directly
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3
  
# output motor_unit config in a readable format
if True:
  import sys
  # only on rank 0
  if (int)(sys.argv[-2]) == 0:
    for mu_no,item in enumerate(motor_units):    
      print("MU {}".format(mu_no))
      for (key,value) in item.items():
        if key != "jitter":
          print("  {}: {}".format(key,value))
  
# boomerAMG literature: https://mooseframework.inl.gov/application_development/hypre.html
# boomeramg options:
#HYPRE preconditioner options
#  -pc_hypre_type <boomeramg> (choose one of) pilut parasails boomeramg ams (PCHYPRESetType)
#HYPRE BoomerAMG Options
#  -pc_hypre_boomeramg_cycle_type <V> (choose one of) V W (None)
#  -pc_hypre_boomeramg_max_levels <25>: Number of levels (of grids) allowed (None)
#  -pc_hypre_boomeramg_max_iter <1>: Maximum iterations used PER hypre call (None)
#  -pc_hypre_boomeramg_tol <0.>: Convergence tolerance PER hypre call (0.0 = use a fixed number of iterations) (None)
#  -pc_hypre_boomeramg_truncfactor <0.>: Truncation factor for interpolation (0=no truncation) (None)
#  -pc_hypre_boomeramg_P_max <0>: Max elements per row for interpolation operator (0=unlimited) (None)
#  -pc_hypre_boomeramg_agg_nl <0>: Number of levels of aggressive coarsening (None)
#  -pc_hypre_boomeramg_agg_num_paths <1>: Number of paths for aggressive coarsening (None)
#  -pc_hypre_boomeramg_strong_threshold <0.25>: Threshold for being strongly connected (None)
#  -pc_hypre_boomeramg_max_row_sum <0.9>: Maximum row sum (None)
#  -pc_hypre_boomeramg_grid_sweeps_all <1>: Number of sweeps for the up and down grid levels (None)
#  -pc_hypre_boomeramg_nodal_coarsen <0>: Use a nodal based coarsening 1-6 (HYPRE_BoomerAMGSetNodal)
#  -pc_hypre_boomeramg_vec_interp_variant <0>: Variant of algorithm 1-3 (HYPRE_BoomerAMGSetInterpVecVariant)
#  -pc_hypre_boomeramg_grid_sweeps_down <1>: Number of sweeps for the down cycles (None)
#  -pc_hypre_boomeramg_grid_sweeps_up <1>: Number of sweeps for the up cycles (None)
#  -pc_hypre_boomeramg_grid_sweeps_coarse <1>: Number of sweeps for the coarse level (None)
#  -pc_hypre_boomeramg_smooth_type <Schwarz-smoothers> (choose one of) Schwarz-smoothers Pilut ParaSails Euclid (None)
#  -pc_hypre_boomeramg_smooth_num_levels <25>: Number of levels on which more complex smoothers are used (None)
#  -pc_hypre_boomeramg_eu_level <0>: Number of levels for ILU(k) in Euclid smoother (None)
#  -pc_hypre_boomeramg_eu_droptolerance <0.>: Drop tolerance for ILU(k) in Euclid smoother (None)
#  -pc_hypre_boomeramg_eu_bj: <FALSE> Use Block Jacobi for ILU in Euclid smoother? (None)
#  -pc_hypre_boomeramg_relax_type_all <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#  -pc_hypre_boomeramg_relax_type_down <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#  -pc_hypre_boomeramg_relax_type_up <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#  -pc_hypre_boomeramg_relax_type_coarse <Gaussian-elimination> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#  -pc_hypre_boomeramg_relax_weight_all <1.>: Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps) (None)
#  -pc_hypre_boomeramg_relax_weight_level <1.>: Set the relaxation weight for a particular level (weight,level) (None)
#  -pc_hypre_boomeramg_outer_relax_weight_all <1.>: Outer relaxation weight for all levels (-k = determined with k CG steps) (None)
#  -pc_hypre_boomeramg_outer_relax_weight_level <1.>: Set the outer relaxation weight for a particular level (weight,level) (None)
#  -pc_hypre_boomeramg_no_CF: <FALSE> Do not use CF-relaxation (None)
#  -pc_hypre_boomeramg_measure_type <local> (choose one of) local global (None)
#  -pc_hypre_boomeramg_coarsen_type <Falgout> (choose one of) CLJP Ruge-Stueben  modifiedRuge-Stueben   Falgout  PMIS  HMIS (None)
#  -pc_hypre_boomeramg_interp_type <classical> (choose one of) classical   direct multipass multipass-wts ext+i ext+i-cc standard standard-wts   FF FF1 (None)
#  -pc_hypre_boomeramg_print_statistics: Print statistics (None)
#  -pc_hypre_boomeramg_print_statistics <3>: Print statistics (None)
#  -pc_hypre_boomeramg_print_debug: Print debug information (None)
#  -pc_hypre_boomeramg_nodal_relaxation: <FALSE> Nodal relaxation via Schwarz (None)
#
