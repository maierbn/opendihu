# This file contains all global variables for the fibers_emg example and their default values. These are the parameters and other internal variables.
# These values will be used by all involved scripts: helper.py, create_partitioned_meshes_for_settings.py and settings_fibers_emg.py
# settings_fibers_emg.py handles setting the parameter values. Those can be overridden on the command line and by specifying a custom variables.py script
# To run the simulation use the settings_fibers_emg.py file, which imports this file, e.g. ./fibers_emg ../settings_fibers_emg.py custom_variables.py

# timing parameters
# -----------------
end_time = 20.0                   # [ms] end time of the simulation
dt_muscle_spindles          = 1e-3  # [ms]

output_timestep_spindles = 1         # [ms] timestep for output of files for all sensor organs and neurons

#--------------------------------
n_elements_muscle = [4,4,10]
meshes = {}
scenario_name = "spindles"

import os
import sys
import numpy as np

input_directory = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input")


# muscle spindles
n_muscle_spindles = 3
muscle_spindle_cellml_file = input_directory+"/mileusenic_muscle_spindle_equilibrium.cellml"
muscle_spindle_mappings = {
  # first: don't define any mappings. then fix all warnings and errors

  # define the order of the parameters for initial values and callbacks
  # opendihu                    cellml
  ("parameter", 0):            "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("parameter", 1):            "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("parameter", 2):            "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("parameter", 3):            "modell/gamma_sta",                 # Static fusimotor drive
  ("parameter", 4):            "modell/gamma_dyn",                 # Dynamic fusimotor drive
  # we first have to define the cellml constants as parameters. primary_afferent ist schon ein state. 
  # opendihu                    cellml
  ("connectorSlot", "ms_out"): "modell/primary_afferent",          # Ia Afferent firing frequency
  ("connectorSlot", "ms_in0"): "modell/L",                         # Spinndle Stretch (= fiber stretch)
  ("connectorSlot", "ms_in1"): "modell/L_dot",                     # Spinndle velocity (= dL/dt)
  ("connectorSlot", "ms_in2"): "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
  ("connectorSlot", "ms_in3"): "modell/gamma_sta",                 # Static fusimotor drive
  ("connectorSlot", "ms_in4"): "modell/gamma_dyn",                 # Dynamic fusimotor drive
}
muscle_spindle_parameters_initial_values = [1, 0, 0, 0, 0]    # [1,0,0,0,0] [L, L_dot, L_ddot, gamma_sta, gamma_dyn] (see above)
muscle_spindle_delay = 30             # [ms] signal delay between muscle spindle model and motoneuron model




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
  # input_values contains a list of λ values from the muscle spindle nodes
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_muscle_spindles (per output slot if there are multiple)
  
  print()
  print()
  print("callback_muscle_spindles_input, input={}".format(input_values), flush=True)
  
  # initialize buffer, buffer is needed to compute velocity and acceleration
  if "stretch" not in buffer:
    buffer["stretch"]      = [1 for _ in range(n_input_values)]
    buffer["velocity"]     = [0 for _ in range(n_input_values)]
    buffer["acceleration"] = [0 for _ in range(n_input_values)]
    buffer["current_time"] = 0
    delta_t = 1
  else:  
    delta_t  = buffer["current_time"] - current_time
    
  # normal stimulation: parse input value of λ
  for i in range(n_input_values):
    stretch = input_values[i]
    
    # the first time when the stretch is not yet computed it has a value of 0, set to 1
    if stretch == 0:
      stretch = 1
      
    # compute velocity and acceleration by difference quotient
    velocity     = (buffer["stretch"][i] - stretch) / delta_t
    acceleration = (buffer["velocity"][i] - velocity) / delta_t
    buffer["stretch"][i]      = stretch
    buffer["velocity"][i]     = velocity
    buffer["acceleration"][i] = acceleration
    print("spindle no. {:2d}, stretch: {:+.8e}, v: {:+.8e}, a: {:+.8e}".format(i, stretch, velocity, acceleration))
      
    # scale value for muscle spindle model 
    output_values[0][i] = stretch                 # fiber stretch, output units: [L0] with L0=cm
    output_values[1][i] = velocity*1e-3           # velocity, output units: [cm*s^-1], 1e-3 cm*ms^-1 = 1e-3 cm*(1e-3s)^-1 = 1e-3 cm*1e3*s^-1 = 1 cm*s^-1
    output_values[2][i] = acceleration*1e-6       # acceleration, output units: [cm*s^-2], 1e-6 cm*ms^-2 = 1e-6 cm*(1e-3s)^-2 = 1e-6 cm*(1e6)*s^-2 = 1 cm*s^-2
    output_values[3][i] = 0                       # static fusimotor drive
    output_values[4][i] = 0                       # dynamic fusimotor drive
  
  # store current time
  buffer["current_time"] = current_time

  sys.stdout.flush() # flush output. neccessary if stdout ist a pipe (eg | tee). files seem to be ok
  
  # artifical muscle spindle input for debugging
  if False:
    for i in range(n_output_values):
      T = 1000   # [ms] cycle duration of stimulus 
      output_values[0][i] = np.sin(current_time / T * np.pi)**2 * 0.001
      output_values[0][i] = 1
      
      # "ms_in0" -> "modell/L",                         # Spinndle Stretch (= fiber stretch)
      # "ms_in1" -> "modell/L_dot",                     # Spinndle velocity (= dL/dt)
      # "ms_in2" -> "modell/L_ddot",                    # Spinndle Acceleration (= d^2L/dt^2)
      # "ms_in3" -> "modell/gamma_sta",                 # Static fusimotor drive
      # "ms_in4" -> "modell/gamma_dyn",                 # Dynamic fusimotor drive
      output_values[0][i] = 1 + np.sin(current_time / T * np.pi) * 0.1                  # Stretch    
      output_values[1][i] = np.cos(current_time / T * np.pi) * 0.1 * (np.pi/T)          # Velocity
      output_values[2][i] = np.sin(current_time / T * np.pi) * (-0.1) * (np.pi/T)**2    # Acceleration
      output_values[3][i] = 0     # static fusimotor drive
      output_values[4][i] = 0     # dynamic fusimotor drive


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
  # mapping muscle spindles output -> motor neuron signals, delay signals by muscle_spindle_delay
  
  # get number of input and output values
  n_input_values = len(input_values)      # = n_muscle_spindles
  n_output_values = len(output_values[0]) # = n_motoneurons

  # initialize buffer the first time, buffer later stores time of last activation, to model signal delay
  if 0 not in buffer:
    for muscle_spindle_index in range(n_input_values):
      buffer[muscle_spindle_index] = None
  
  # total summed up signal
  total_signal = 0
      
  # loop over muscle spindles
  for muscle_spindle_index in range(n_input_values):
    
    # determine spike by threshold
    if input_values[muscle_spindle_index] > 0:
      buffer[muscle_spindle_index] = current_time    # store time of last activation in buffer
      
    # if there has been a stimulation so far
    if buffer[muscle_spindle_index] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = muscle_spindle_delay             # [ms] delay of the signal
      gaussian_std_dev = 10                     # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[muscle_spindle_index]) * 5
        
      # sum up all input signals
      total_signal += delayed_signal * 1e-3
    
      print(" spindle {:2d}/{:2d} signal: {:.2e}".format(muscle_spindle_index,n_input_values,delayed_signal * 1e-3))
        
  # compute average signal over all muscle spindles
  total_signal /= n_input_values
    
  # motor neuron fires with ~14Hz if drive(t) = 5e-3
  
  # set values to all connected motoneurons
  for motoneuron_index in range(n_output_values):
    output_values[0][motoneuron_index] = total_signal
  
  print("muscle_spindles_to_motoneurons: {} -> {}".format(input_values, output_values))
  sys.stdout.flush() # flush output. neccessary if stdout ist a pipe (eg | tee). files seem to be ok
