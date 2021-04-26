#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script to visualize neuron firings.
#

import sys, os
import numpy as np
import time
import pickle
import py_reader    # reader utility for opendihu *.py files

show_plot = True

# import needed packages from matplotlib
import matplotlib as mpl
if not show_plot:
  mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec

# get all input data in current directory
filenames = os.listdir("out")

# collect the filenames
condition = lambda filename: ".py" in filename and "muscle_spindles" in filename
muscle_spindles_files = ["out/{}".format(filename) for filename in sorted(list(np.extract(np.array(list(map(condition, filenames))), filenames)))]

condition = lambda filename: ".py" in filename and "golgi_tendon_organs" in filename
golgi_tendon_organs_files = ["out/{}".format(filename) for filename in sorted(list(np.extract(np.array(list(map(condition, filenames))), filenames)))]

condition = lambda filename: ".py" in filename and "interneurons" in filename
interneurons_files = ["out/{}".format(filename) for filename in sorted(list(np.extract(np.array(list(map(condition, filenames))), filenames)))]

condition = lambda filename: ".py" in filename and "motoneurons" in filename
motoneurons_files = ["out/{}".format(filename) for filename in sorted(list(np.extract(np.array(list(map(condition, filenames))), filenames)))]

print("Number of input files: muscle spindles: {}, golgi tendon organs: {}, interneurons: {}, motoneurons: {}".\
  format(len(muscle_spindles_files), len(golgi_tendon_organs_files), len(interneurons_files), len(motoneurons_files)))


# load data
muscle_spindles_data = py_reader.load_data(muscle_spindles_files)
golgi_tendon_organs_data = py_reader.load_data(golgi_tendon_organs_files)
interneurons_data = py_reader.load_data(interneurons_files)
motoneurons_data = py_reader.load_data(motoneurons_files)


# create plots
fig,axes = plt.subplots(4,1,figsize=(12,8),sharex=True)

# ---------------------
# plot muscle spindles
component_name_input = "(P)membrane/i_Stim"
component_name_output = "membrane/V"

t_values = None
values_output = None
values_input = None

# loop over datasets at different times
for i,dataset in enumerate(muscle_spindles_data):
  
  # get the data for the current timestep
  data_input = py_reader.get_values(dataset, "parameters", component_name_input)
  data_output = py_reader.get_values(dataset, "solution", component_name_output)
  
  if data_input is None:
    print("No data found for Golgi tendon organs or component '{}' does not exist.\n".format(component_name_input))
  if data_output is None:
    print("No data found for Golgi tendon organs or component '{}' does not exist.\n".format(component_name_output))
  
  # create arrays the first time
  if values_output is None:
    values_input = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    values_output = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    t_values = np.zeros((len(golgi_tendon_organs_data)))
  
  # store values
  values_input[:,i] = data_input
  values_output[:,i] = data_output
  t_values[i] = dataset['currentTime']

# plot lines for all timesteps
# loop over neurons
n = values_output.shape[0]
for i in range(n):
  color = next(axes[0]._get_lines.prop_cycler)['color']
  axes[0].plot(t_values, values_output[i,:], '-', color=color)
  if i == 0:
    ax2 = axes[0].twinx()
  ax2.plot(t_values, values_input[i,:], ':', color=color)

# set title and axis labels
axes[0].set_title('Muscle spindles (number: {})'.format(n))
axes[0].set_ylabel('voltage [mV]\n(solid lines)')
ax2.set_ylabel('input current [uA]\n(dotted lines)')

# ---------------------
# plot Golgi tendon organs
component_name_input = "(P)membrane/i_Stim"
component_name_output = "membrane/V"

t_values = None
values_input = None
values_output = None

# loop over datasets at different times
for i,dataset in enumerate(golgi_tendon_organs_data):
  
  # get the data for the current timestep
  data_input = py_reader.get_values(dataset, "parameters", component_name_input)
  data_output = py_reader.get_values(dataset, "solution", component_name_output)
  
  if data_input is None:
    print("No data found for Golgi tendon organs or component '{}' does not exist.\n".format(component_name_input))
  if data_output is None:
    print("No data found for Golgi tendon organs or component '{}' does not exist.\n".format(component_name_output))
  
  # create arrays the first time
  if values_output is None:
    values_input = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    values_output = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    t_values = np.zeros((len(golgi_tendon_organs_data)))
  
  # store values
  values_input[:,i] = data_input
  values_output[:,i] = data_output
  t_values[i] = dataset['currentTime']

# plot lines for all timesteps
# loop over neurons
n = values_output.shape[0]
for i in range(n):
  color = next(axes[1]._get_lines.prop_cycler)['color']
  axes[1].plot(t_values, values_output[i,:], '-', color=color)
  if i == 0:
    ax2 = axes[1].twinx()
  ax2.plot(t_values, values_input[i,:], ':', color=color)

# set title and axis labels
axes[1].set_title('Golgi tendon organs (number: {})'.format(n))
axes[1].set_ylabel('voltage [mV]\n(solid lines)')
ax2.set_ylabel('input current [uA]\n(dotted lines)')

# ---------------------
# plot interneurons
component_name_input = "(P)membrane/i_Stim"
component_name_output = "membrane/V"

t_values = None
values_output = None
values_input = None

# loop over datasets at different times
for i,dataset in enumerate(interneurons_data):
  
  # get the data for the current timestep
  data_input = py_reader.get_values(dataset, "parameters", component_name_input)
  data_output = py_reader.get_values(dataset, "solution", component_name_output)
  
  if data_input is None:
    print("No data found for interneurons or component '{}' does not exist.\n".format(component_name_input))
  if data_output is None:
    print("No data found for interneurons or component '{}' does not exist.\n".format(component_name_output))
  
  # create arrays the first time
  if values_output is None:
    values_input = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    values_output = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    t_values = np.zeros((len(golgi_tendon_organs_data)))
  
  # store values
  values_input[:,i] = data_input
  values_output[:,i] = data_output
  t_values[i] = dataset['currentTime']

# plot lines for all timesteps
# loop over neurons
n = values_output.shape[0]
for i in range(n):
  color = next(axes[2]._get_lines.prop_cycler)['color']
  axes[2].plot(t_values, values_output[i,:], '-', color=color)
  if i == 0:
    ax2 = axes[2].twinx()
  ax2.plot(t_values, values_input[i,:], ':', color=color)
  

# set title and axis labels
axes[2].set_title('Interneurons (number: {})'.format(n))
axes[2].set_ylabel('voltage [mV]\n(solid lines)')
ax2.set_ylabel('input current [uA]\n(dotted lines)')


# ---------------------
# plot motoneurons
component_name_input = "(P)motor_neuron/drive"
component_name_output = "motor_neuron/V_s"

t_values = None
values_output = None
values_input = None

# loop over datasets at different times
for i,dataset in enumerate(motoneurons_data):
  
  # get the data for the current timestep
  data_input = py_reader.get_values(dataset, "parameters", component_name_input)
  data_output = py_reader.get_values(dataset, "solution", component_name_output)
  
  if data_input is None:
    print("No data found for motoneurons or component '{}' does not exist.\n".format(component_name_input))
  if data_output is None:
    print("No data found for motoneurons or component '{}' does not exist.\n".format(component_name_output))
  
  # create arrays the first time
  if values_output is None:
    values_input = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    values_output = np.zeros((len(data_output), len(golgi_tendon_organs_data)))   # each column is the data for one timestep, for multiple neurons
    t_values = np.zeros((len(golgi_tendon_organs_data)))
  
  # store values
  values_input[:,i] = data_input
  values_output[:,i] = data_output
  t_values[i] = dataset['currentTime']

# plot lines for all timesteps
# loop over neurons
n = values_output.shape[0]
for i in range(n):
  color = next(axes[3]._get_lines.prop_cycler)['color']
  axes[3].plot(t_values, values_output[i,:], '-', color=color)
  if i == 0:
    ax2 = axes[3].twinx()
  ax2.plot(t_values, values_input[i,:], ':', color=color)
  
# set title and axis labels
axes[3].set_title('Motor neurons (number: {})'.format(n))
axes[3].set_xlabel('time [ms]')
axes[3].set_ylabel('voltage [mV]\n(solid lines)')
ax2.set_ylabel('input current [uA]\n(dotted lines)')

# show plot window
plt.savefig("plot.png")
plt.show()
