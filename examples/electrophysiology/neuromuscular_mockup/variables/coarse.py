
# scenario name for log file
scenario_name = "coarse"

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
bottom_traction = [0.0,0.0,0.0]        # [N]

# timing and activation parameters
# -----------------
import random
random.seed(0)  # ensure that random numbers are the same on every rank

n_compartments = 1

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
snes_max_iterations = 10                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 2       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
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
end_time = 4000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver, i.e. the diffusion
dt_splitting = dt_multidomain       # [ms] timestep width of strang splitting between 0D and multidomain, this is the same as the dt_multidomain, because we do not want to subcycle for the diffusion part
dt_elasticity = 1e-1                # [ms] time step width of elasticity solver

dt_neurons = 1e-3                   # [ms] time step width for all neuron solvers
dt_muscle_spindles = 1e-3           # [ms] timestep width of cellml solver of muscle spindles
dt_golgi_tendon_organs = 1e-3       # [ms] timestep width of cellml solver of golgi tendon organs
dt_muscle_spindles = 1e-3           # [ms] timestep width of cellml solver of muscle spindles
dt_interneuron = 1e-3               # [ms] timstep width of the cellml solver for interneurons
dt_motoneuron = 1e-3                # [ms] timstep width of the cellml solver for motoneurons
dt_stimulation_check = 1e-3         # [ms] interval when to apply stimulation input to electro-mechanics solver

output_timestep_multidomain = 1e-1  # [ms] timestep for fiber output, 0.5
output_timestep_elasticity = 1      # [ms] timestep for elasticity output files
output_timestep_neurons = 1e-3      # [ms] timestep for output of files for neurons
output_timestep_motoneuron = 1e-3   # [ms] timestep for output of files for motoneuron

# input files
#multidomain_cellml_file = "../../input/hodgkin_huxley-razumova.cellml"
fiber_file = "../../input/left_biceps_brachii_9x9fibers.bin"
fiber_file = "../../input/left_biceps_brachii_13x13fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
firing_times_file = "../../input/MU_firing_times_once.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "../../input/MU_fibre_distribution_10MUs.txt"

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
# If you change this, delete the compartment_relative_factors.* files, they have to be generated again.
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 50
sampling_stride_fat = 1

# how much of the multidomain mesh is used for elasticity
sampling_factor_elasticity_x = 0.5    
sampling_factor_elasticity_y = 0.5
sampling_factor_elasticity_z = 0.3
sampling_factor_elasticity_fat_y = 0.5

# neurons and sensors
# muscle spindles
n_muscle_spindles = 6
muscle_spindle_cellml_file = "../../input/hodgkin_huxley_1952.cellml"
muscle_spindle_mappings = {
  ("parameter", 0):           "membrane/i_Stim",   # stimulation
  ("outputConnectorSlot", 0): "membrane/V",        # voltage
  ("outputConnectorSlot", 1): "membrane/i_Stim",   # stimulation
}
muscle_spindle_parameters_initial_values = [0]    # [i_Stim]

# golgi tendon organs
n_golgi_tendon_organs = 4
golgi_tendon_organ_cellml_file = "../../input/hodgkin_huxley_1952.cellml"
golgi_tendon_organ_mappings = {
  ("parameter", 0):           "membrane/i_Stim",   # stimulation
  ("outputConnectorSlot", 0): "membrane/V",        # voltage
  ("outputConnectorSlot", 1): "membrane/i_Stim",   # stimulation
}
golgi_tendon_organ_parameters_initial_values = [0]    # [i_Stim]

# inter neurons
n_interneurons = 6
interneuron_cellml_file = "../../input/hodgkin_huxley_1952.cellml"
interneuron_mappings = {
  ("parameter", 0):           "membrane/i_Stim",   # stimulation
  ("outputConnectorSlot", 0): "membrane/V",        # voltage
  ("outputConnectorSlot", 1): "membrane/i_Stim",   # stimulation
}
interneuron_parameters_initial_values = [0]    # [i_Stim]

# motor neurons
n_motoneurons = 4
motoneuron_cellml_file = "../../input/WSBM_1457_MN_Cisi_Kohn_2008.cellml"
motoneuron_mappings = {
  ("parameter", 0):           "motor_neuron/drive",   # stimulation
  ("outputConnectorSlot", 0): "motor_neuron/V_s",     # voltage
}
motoneuron_parameters_initial_values = [0.01]    # [drive]

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
states_output = True    # if also the subcellular states should be output, this produces large files, set output_timestep_0D_states
show_linear_solver_output = True    # if every solve of multidomain diffusion should be printed
disable_firing_output = True   # if information about firing of MUs should be printed

# callbacks

def callback_muscle_spindles_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)

def callback_muscle_spindles_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)
  
  # initialize buffer the first time
  if 0 not in buffer:
    buffer[0] = None
  
  # determine spike by threshold
  if input_values[0] > 20:
    buffer[0] = current_time    # store time of last activation in buffer[0]
    
  # if there has been a stimulation so far
  if buffer[0] is not None:
    
    # convolute Dirac delta, kernel is a shifted and scaled gaussian
    t_delay = 10              # [ms] delay of the signal
    gaussian_std_dev = 0.1    # [ms] width of the gaussian curve
    convolution_kernel = lambda t: scipy.stats.norm.pdf(t, loc=t_delay, scale=gaussian_std_dev)*np.sqrt(2*np.pi)*gaussian_std_dev
    delayed_signal = convolution_kernel(current_time - buffer[0]) * 20
      
    # loop over output values and set all to the computed signal, cut off at 1e-5
    if delayed_signal > 1e-5:
      print("motoneuron t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[0], delayed_signal))
      for i in range(n_output_values):
        output_values[i] = delayed_signal
    else:
      for i in range(n_output_values):
        output_values[i] = None     # do not set any values
    

def callback_golgi_tendon_organs_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)

def callback_golgi_tendon_organs_to_interneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)

def callback_interneurons_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)

def callback_motoneurons_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)
  
def callback_motoneuron_output(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of float values) Initially, this is a list of the form [None, None, ..., None] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # get number of input and output values
  n_input_values = len(input_values)      # =1 here (1 motoneuron)
  n_output_values = len(output_values)    # =3 here (3 points in neuromuscular junction)
  
