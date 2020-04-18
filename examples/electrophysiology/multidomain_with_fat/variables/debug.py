# debugging scenario that produces a fast simulation, not really physiological conditions (Am is too high therefore more diffusion than normal)
# scenario name for log file
scenario_name = "debug"

# material parameters
# --------------------
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
sigma_xf = 0                # [mS/cm] conductivity in cross-fiber direction (xf)
sigma_e_f = 6.7             # [mS/cm] conductivity in extracellular space, fiber direction (f)
sigma_e_xf = 3.35           # [mS/cm] conductivity in extracellular space, cross-fiber direction (xf) / transverse

Am = 500.0                  # [cm^-1] surface area to volume ratio, actual values will be set by motor_units, not by this variable
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch), actual values will be set by motor_units, not by this variable

Am = 1.0  # set Am lower such that there is more diffusion
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor unit parameters
# stimulation frequency in [Hz], activation_start_time in [s]
#   fiber_no: center MU around this fiber, create f_r as gaussion from standard_deviation and maximum

motor_units = [
  {"fiber_no": 10, "standard_deviation": 0.2, "maximum": 0.5, "activation_start_time": 0.0, "stimulation_frequency": 10.0, "cm": 0.58},
  {"fiber_no": 30, "standard_deviation": 0.2, "maximum": 0.4, "activation_start_time": 0.0, "stimulation_frequency": 10.0, "cm": 0.58},
  {"fiber_no": 40, "standard_deviation": 0.3, "maximum": 0.6, "activation_start_time": 0.0, "stimulation_frequency": 10.0, "cm": 0.58},
]
motor_units = motor_units[0:2]  # only two motor units

# solvers
# -------
multidomain_solver_type = "cg"         # solver for the multidomain problem
multidomain_preconditioner_type = "euclid"   # preconditioner
multidomain_alternative_solver_type = "lu"           # alternative solver, used when normal solver diverges
multidomain_alternative_preconditioner_type = "none"     # preconditioner of the alternative solver

theta = 0.5                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True    # if the diagonal blocks of the system matrix should be used as preconditioner matrix, (True means more iterations in the solver but faster preconditioner)
use_lumped_mass_matrix = False             # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations

# set initial guess to zero for direct solver
initial_guess_nonzero = "lu" not in multidomain_solver_type 

multidomain_absolute_tolerance = 1e-15    # absolute residual tolerance for the multidomain solver
multidomain_relative_tolerance = 1e-15    # absolute residual tolerance for the multidomain solver
theta = 0.5                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True    # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations

# timing parameters
# -----------------
end_time = 4000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs (3e-3)
dt_multidomain = dt_0D              # [ms] timestep width of multidomain
dt_splitting = dt_multidomain       # [ms] overall timestep width of strang splitting (3e-3)
output_timestep_multidomain = 1e-1  # [ms] timestep for multidomain output

# input files
cellml_file = "../../input/hodgkin_huxley_1952.c"
fiber_file = "../../input/left_biceps_brachii_7x7fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../input/MU_firing_times_immediately.txt"


# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 74   # faster, but stimulus does not propagate
sampling_stride_fat = 2


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
  return motor_units[mu_no % len(motor_units)]["cm"]
  #return Cm
  
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
