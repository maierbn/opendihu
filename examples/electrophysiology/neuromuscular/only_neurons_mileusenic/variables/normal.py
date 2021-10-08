
# scenario name for log file
scenario_name = "normal"

# timing and activation parameters
# -----------------
import random
random.seed(0)  # ensure that random numbers are the same on every rank
import numpy as np

# timing parameters
# -----------------
end_time = 5_000.0                  # [ms] end time of the simulation
#end_time = 100.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of the multidomain solver, i.e. the diffusion
dt_splitting = dt_multidomain       # [ms] timestep width of strang splitting between 0D and multidomain, this is the same as the dt_multidomain, because we do not want to subcycle for the diffusion part
dt_elasticity = 1                   # [ms] timestep width of elasticity solver

dt_neurons = 1e-2                   # [ms] same timestep width for all neuron solvers
dt_golgi_tendon_organs = dt_neurons # [ms] timestep width of cellml solver of golgi tendon organs
dt_muscle_spindles     = dt_neurons # [ms] timestep width of cellml solver of muscle spindles
dt_interneuron         = dt_neurons # [ms] timestep width of the cellml solver for interneurons
dt_motoneuron          = dt_neurons # [ms] timestep width of the cellml solver for motoneurons

dt_neuron_transfer     = 10*dt_neurons  # [ms] interval when to call callback functions and transfer values between CellML models, increase this to speed up the simulation

output_timestep_multidomain = 5     # [ms] timestep for multidomain solver output
output_timestep_elasticity = 5      # [ms] timestep for elasticity output files
output_timestep_neurons = 5         # [ms] timestep for output of files for all sensor organs and neurons
output_timestep_motoneuron = 5      # [ms] timestep for output of files for motoneuron

# input files
# -----------
# neurons and sensors
# muscle spindles
n_muscle_spindles = 1
muscle_spindle_cellml_file = "../../../input/hodgkin_huxley_1952.cellml"
muscle_spindle_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("connectorSlot", "ms_out"): "membrane/V",        # voltage
  ("connectorSlot", "ms_in"):  "membrane/i_Stim",   # stimulation
}
muscle_spindle_parameters_initial_values = [0]    # [i_Stim]
muscle_spindle_delay = 30             # [ms] signal delay between muscle spindle model and motoneuron model

# golgi tendon organs
n_golgi_tendon_organs = 3
golgi_tendon_organ_cellml_file = "../../../input/hodgkin_huxley_1952.cellml"
golgi_tendon_organ_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("connectorSlot", "gt_out"): "membrane/V",        # voltage
  ("connectorSlot", "gt_in"):  "membrane/i_Stim",   # stimulation
}
golgi_tendon_organ_parameters_initial_values = [0]    # [i_Stim]
golgi_tendon_organ_delay = 300

# interneurons
n_interneurons = 3
interneuron_cellml_file = "../../../input/hodgkin_huxley_1952.cellml"
interneuron_mappings = {
  ("parameter", 0):            "membrane/i_Stim",   # stimulation
  ("parameter", 1):            "membrane/Cm",       # stimulation
  ("connectorSlot", "in_out"): "membrane/V",        # voltage
  ("connectorSlot", "in_in"):  "membrane/i_Stim",   # stimulation
}
# initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
interneuron_parameters_initial_values = []    # [i_Stim, Cm]
for i in range(n_interneurons):
  
  # compute a scaling factor that runs exponentially from min_factor to max_factor
  min_factor = 0.5
  max_factor = 2.0
  
  # ansatz scaling_factor(i) = c1 + c2*exp(i),
  # scaling_factor(0) = min = c1 + c2  =>  c1 = min - c2
  # scaling_factor(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_factor - min_factor) / (np.exp(n_interneurons-1) - 1)
  c1 = min_factor - c2
  scaling_factor = c1 + c2*np.exp(i)
  
  # add parameter values for motoneuron i
  interneuron_parameters_initial_values += [0.0, 1*scaling_factor]
  
print(interneuron_parameters_initial_values)
  
# motor neurons
n_motoneurons = 3
motoneuron_cellml_file = "../../../input/WSBM_1457_MN_Cisi_Kohn_2008.cellml"
motoneuron_mappings = {
  ("parameter", 0):            "motor_neuron/drive",   # stimulation
  ("parameter", 1):            "lumped_geometry_parameters/C_m",
  ("parameter", 2):            "lumped_geometry_parameters/R_md",
  ("connectorSlot", "mn_out"): "motor_neuron/V_s",     # voltage
  ("connectorSlot", "mn_in"):  "motor_neuron/drive",   # stimulation
}
# initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
motoneuron_parameters_initial_values = []
for i in range(n_motoneurons):
  
  # compute a scaling factor that runs exponentially from min_factor to max_factor
  min_factor = 0.9
  max_factor = 2.0
  
  # ansatz scaling_factor(i) = c1 + c2*exp(i),
  # scaling_factor(0) = min = c1 + c2  =>  c1 = min - c2
  # scaling_factor(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_factor - min_factor) / (np.exp(n_motoneurons-1) - 1)
  c1 = min_factor - c2
  scaling_factor = c1 + c2*np.exp(i)
  
  # add parameter values for motoneuron i
  motoneuron_parameters_initial_values += [0.01, 1*scaling_factor, 14.4*scaling_factor]

# other options
n_compartments = 1

# callbacks

def callback_muscle_spindles_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # map from λ in the 3D mesh to muscle spindles model input
  
  # get number of input and output values
  n_input_values = len(input_values)      # = 0 (normally n_muscle_spindles, if the muscle mesh is connected)
  n_output_values = len(output_values[0]) # = n_muscle_spindles (per output slot if there are multiple)
  
  # normal stimulation, will not be executed here, because n_input_values is 0
  for i in range(n_input_values):
    stretch = input_values[i]
    
    # the first time when the stretch is not yet computed it has a value of 0, set to 1
    if stretch == 0:
      stretch = 1
      
    output_values[0][i] = abs(stretch-1) * 150
  
  # artifical muscle spindle input
  for i in range(n_output_values):
    T = 1000   # [ms] cycle duration of stimulus 
    output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 10
  
  #print("stretch at muscle spindles: {}, output: {}".format(input_values, output_values))
  
def callback_muscle_spindles_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping muscle spindles output -> motor neuron signals, delay signals by 
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_muscle_spindles (per output slot if there are multiple)
  
  # initialize buffer the first time
  if 0 not in buffer:
    for muscle_spindle_index in range(n_input_values):
      buffer[muscle_spindle_index] = None
  
  # loop over muscle spindles
  for muscle_spindle_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[muscle_spindle_index] > 0:
      buffer[muscle_spindle_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[muscle_spindle_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = muscle_spindle_delay              # [ms] delay of the signal
      gaussian_std_dev = 10                      # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[muscle_spindle_index]) * 5
        
      # loop over output values and set all to the computed signal, cut off at 1e-5
      if delayed_signal > 1e-5:
        #print("muscle spindle t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[muscle_spindle_index], delayed_signal))
        output_values[0][muscle_spindle_index] = delayed_signal
      else:
        output_values[0][muscle_spindle_index] = None     # do not set any values
        
  #print("muscle_spindles_to_motoneurons: {} -> {}".format(input_values, output_values))

def callback_golgi_tendon_organs_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # map from λ in the 3D mesh to golgi tendon organs
  
  # get number of input and output values
  n_input_values = len(input_values)      # = 0 (normally n_golgi_tendon_organs, if muscle mesh is connected)
  n_output_values = len(output_values[0]) # = n_golgi_tendon_organs (per output slot if there are multiple)
  
  # normal stimulation, will not be executed here, because n_input_values is 0
  for i in range(n_input_values):
    stretch = input_values[i]
    
    # the first time when the stretch is not yet computed it has a value of 0, set to 1
    if stretch == 0:
      stretch = 1
      
    output_values[0][i] = abs(stretch-1) * 200
    
  # artifical muscle spindle input
  for i in range(n_output_values):
    T = 100   # [ms] cycle duration of stimulus 
    
    if 2000 + i*200 < current_time  <= 2100 + i*200:
      output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 20
  
  #print("stretch at Golgi tendon organs: {}, output: {}".format(input_values, output_values))

def callback_golgi_tendon_organs_to_interneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping Golgi tendon organs -> interneurons
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_golgi_tendon_organs
  n_output_values = len(output_values[0]) # = n_interneurons (per output slot if there are multiple)
  
  # sum up all input signals and set all outputs to this sum
  
  # collect sum of all Golgi tendon organs
  total_signal = 0
  
  # loop over Golgi tendon organs
  for golgi_tendon_organ_index in range(n_input_values):
    
    # if input is active
    if input_values[golgi_tendon_organ_index] > 20:
      total_signal += input_values[golgi_tendon_organ_index]
    
  # set same value to all connected interneurons
  for interneuron_index in range(n_output_values):
    output_values[0][interneuron_index] = total_signal * 0.01
    
  # initialize buffer the first time
  if 0 not in buffer:
    for golgi_tendon_organ_index in range(n_input_values):
      buffer[golgi_tendon_organ_index] = None
  
  # loop over golgi tendon organ inputs
  for golgi_tendon_organ_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[golgi_tendon_organ_index] > 20:
      buffer[golgi_tendon_organ_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[golgi_tendon_organ_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = 0                   # [ms] delay of the signal
      gaussian_std_dev = 10         # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[golgi_tendon_organ_index]) * 5
        
      output_values[0][golgi_tendon_organ_index] = delayed_signal
  
  #print("golgi_tendon_organs_to_interneurons input: {}, output: {}".format(input_values, output_values))

def callback_interneurons_to_motoneurons(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  # mapping interneurons -> input for motor neurons, i.e. signal delay from interneurons to motoneurons
  # the actual N->M mapping from interneurons to motoneurons is done by callback_motoneurons_input
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_interneurons
  n_output_values = len(output_values[0]) # = n_interneurons (per output slot if there are multiple)
  
  # initialize buffer the first time
  if 0 not in buffer:
    for interneuron_index in range(n_input_values):
      buffer[interneuron_index] = None
  
  # loop over interneurons
  for interneuron_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[interneuron_index] > 0:
      buffer[interneuron_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[interneuron_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = golgi_tendon_organ_delay          # [ms] delay of the signal
      gaussian_std_dev = 10                       # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[interneuron_index]) * 5
        
      # loop over output values and set all to the computed signal, cut off at 1e-5
      if delayed_signal > 1e-5:
        #print("interneuron t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[interneuron_index], delayed_signal))
        output_values[0][interneuron_index] = delayed_signal
      else:
        output_values[0][interneuron_index] = None     # do not set any values
  
  #print("interneurons_to_motoneurons input: {}, output: {}".format(input_values, output_values))
  
def callback_motoneurons_input(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # map from delayed muscle spindle model outputs and delayed interneuron outputs to motoneuron inputs
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles + n_interneurons
  n_output_values = len(output_values[0]) # = n_motoneurons (per output slot if there are multiple)
  
  # sum up all input signals and set all outputs to this sum
  
  # collect sum of all Golgi tendon organs
  total_signal = 0
  
  # loop over input values
  for input_index in range(n_input_values):
    total_signal += input_values[input_index] * 1e-3
    
  # add cortical input
  total_signal += 2e-3
  
  # set same value to all connected motoneurons
  for motoneuron_index in range(n_output_values):
    output_values[0][motoneuron_index] = total_signal
    
  #print("motoneurons input: {}, output: {}".format(input_values, output_values))
  
def callback_motoneuron_output(input_values, output_values, current_time, slot_nos, buffer):
  """
  Callback function that transform a number of input_values to a number of output_values.
  This function gets called by a MapDofs object.
  :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
  :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                        The output values buffer, potentially for multiple slots.
                        Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                        the number of required output values. The function should set some of the entries to a computed value.
                        The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
  :param current_time:  Current simulation time.
  :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
  :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                        this function every time. Using this buffer, it is possible to implement a time delay of signals.
  """
  
  # This is callback is currently not used, because the MapDofs object has the mode set to "localSetIfAboveThreshold",
  # which assigns the value 20 to all outputs if the input is above the threshold of 20
  
  # get number of input and output values
  n_input_values = len(input_values)      # 1 (1 motoneuron)
  n_output_values = len(output_values[0]) # =N (number of points in neuromuscular junction, this is all nodes in the x-y plane at the center of the muscle)
  