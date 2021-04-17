
# scenario name for log file
scenario_name = "coarse"

# Fixed units in cellMl models:
# These define the unit system.
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 uA = 1e-6 A
# 1 uF = 1e-6 F
# 
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
bottom_traction = [0.0,0.0,-1.0]        # [N]

# Monodomain parameters
# --------------------
# quantities in CellML unit system
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
sigma_xf = 0                # [mS/cm] conductivity in cross-fiber direction (xf)
sigma_e_f = 6.7             # [mS/cm] conductivity in extracellular space, fiber direction (f)
sigma_e_xf = 3.35           # [mS/cm] conductivity in extracellular space, cross-fiber direction (xf) / transverse

Am = 500.0                  # [cm^-1] surface area to volume ratio, actual values will be set by motor_units, not by this variable
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch), actual values will be set by motor_units, not by this variable

# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor unit parameters similar to paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
# however, values from paper fail for mu >= 6, then stimulus gets reflected at the ends of the muscle, therefore fiber radius is set to <= 55

import random
random.seed(0)  # ensure that random numbers are the same on every rank
import numpy as np

n_fibers_in_fiber_file = 169
n_motor_units = 100   # number of motor units

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

  # stimulation current for motoneuron (reverse exponential distribution)
  min_value = 2
  max_value = 10
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  mn_istim = c1 + c2*1.02**(n_motor_units-1-mu_no)

  # exponential distribution: low number of fibers per MU, slow twitch (type I), activated first --> high number of fibers per MU, fast twitch (type II), activated last
  motor_units.append(
  {
    "fiber_no":              random.randint(0,n_fibers_in_fiber_file),  # [-] fiber from input files that is the center of the motor unit domain
    "maximum":               maximum,                # [-] maximum value of f_r, create f_r as gaussian with standard_deviation and maximum around the fiber g
    "standard_deviation":    standard_deviation,     # [-] standard deviation of f_r
    "radius":                radius,                 # [μm] parameter for motor unit: radius of the fiber, used to compute Am
    "cm":                    cm,                     # [uF/cm^2] parameter Cm
    "mn_istim":              mn_istim,               # [uA] constant stimulation current for the motoneuron
  })


#motor_units = motor_units[0:3]  # only use 3 motor units
motoneuron_mappings = {
  ("parameter", 0):           "membrane/i_Stim",  
  ("connectorSlot", "vm_mn"): "membrane/V",
}

# set values for parameters: [i_Stim, V_threshold, V_firing], 
# i_Stim = stimulation current of motoneuron, V_threshold = threshold if trans-membrane voltage is above, motoneuron fires, V_firing = value of Vm to set if motoneuron fires
motoneuron_parameters_initial_values = [motor_units[mu_no]["mn_istim"] for mu_no in range(n_motor_units)]

# solvers
# -------
# potential flow
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "none" # preconditioner

# multidomain
multidomain_solver_type = "gmres"          # solver for the multidomain problem
multidomain_preconditioner_type = "euclid"   # preconditioner

multidomain_alternative_solver_type = "gmres"            # alternative solver, used when normal solver diverges
multidomain_alternative_preconditioner_type = "euclid"    # preconditioner of the alternative solver

multidomain_absolute_tolerance = 1e-15 # absolute residual tolerance for the multidomain solver
multidomain_relative_tolerance = 1e-15 # absolute residual tolerance for the multidomain solver

# elasticity
elasticity_solver_type = "lu"
elasticity_preconditioner_type = "none"
snes_max_iterations = 34                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 5       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-5      # absolute tolerance of the nonlinear solver
relative_tolerance = 1e-5           # relative tolerance of the residual of the linear solver
absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver

# set initial guess to zero for direct solver
initial_guess_nonzero = "lu" not in multidomain_solver_type 

theta = 1.0                               # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
use_symmetric_preconditioner_matrix = True   # if the diagonal blocks of the system matrix should be used as preconditioner matrix
use_lumped_mass_matrix = False            # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations

# timing parameters
# -----------------
end_time = 3e-2                   # [ms] end time of the simulation
end_time = 4000                       # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_motoneuron = 1e-3                # [ms] timestep width for motoneurons 
dt_stimulation_check = 1e-2         # [ms] timestep width for when to check if the motoneurons stimulated 
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver, i.e. the diffusion
dt_splitting = dt_multidomain       # [ms] timestep width of strang splitting between 0D and multidomain, this is the same as the dt_multidomain, because we do not want to subcycle for the diffusion part
dt_elasticity = 1e-2                # [ms] time step width of elasticity solver
output_timestep_multidomain = 2 #2e0  # [ms] timestep for fiber output, 0.5
output_timestep_motoneuron = 1e0   # [ms] timestep for output of motoneuron
output_timestep_elasticity = 2     # [ms] timestep for elasticity output files

# input files
#motoneuron_cellml_file = "../../../input/motoneuron_hodgkin_huxley.cellml"
motoneuron_cellml_file = "../../../input/hodgkin_huxley_1952.cellml"
#cellml_file = "../../../input/new_slow_TK_2014_12_08.c"
cellml_file = "../../../input/hodgkin_huxley-razumova.cellml"

#fiber_file = "../../../input/left_biceps_brachii_9x9fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_13x13fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
#firing_times_file = "../../../input/MU_firing_times_once.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "../../../input/MU_fibre_distribution_10MUs.txt"

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
# If you change this, delete the compartment_relative_factors.* files, they have to be generated again.
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 22      # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40
sampling_stride_fat = 1

# how much of the multidomain mesh is used for elasticity
sampling_factor_elasticity_x = 0.5    
sampling_factor_elasticity_y = 0.5
sampling_factor_elasticity_z = 0.3
sampling_factor_elasticity_fat_y = 0.5

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
  
# output motor_unit config in a readable format
if True:
  import sys
  # only on rank 0
  if (int)(sys.argv[-2]) == 0:
    import pprint 
    pp = pprint.PrettyPrinter()
    pp.pprint(motor_units)
