#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to visualize neuron firings.
#

import sys, os
import numpy as np
import time
import pickle
import argparse
import py_reader    # reader utility for opendihu *.py files

show_plot = True

# import needed packages from matplotlib
import matplotlib as mpl
if not show_plot:
  mpl.use('Agg')

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec

# get all input data in current directory
filenames = os.listdir("out")

parallel = False
suffix = ".py"
if parallel:
  suffix = ".0.py"

parser = argparse.ArgumentParser(description='plot_spindles')
parser.add_argument('--store_frequency_file', default="")
parser.add_argument('--load_frequency_file', default="")
parser.add_argument('--plot_firing_times', default=False, action="store_true")
parser.add_argument('--plot_frequency_progression', default=False, action="store_true")
parser.add_argument('--end_time', type=float, default=np.inf)
parser.add_argument('--start_time', type=float, default=0)
parser.add_argument('--files_stride', type=int, default=1)
args = parser.parse_args()

if args.plot_firing_times:
  print("plot firing times")
  
if args.plot_frequency_progression:
  print("plot frequency progression")

if args.end_time < np.inf or args.start_time > 0:
  print("restrict data to [{},{}]".format(args.start_time,args.end_time))

print("store_frequency_file: \"{}\"\nload_frequency_file: \"{}\"".format(args.store_frequency_file, args.load_frequency_file))

# collect the filenames
condition = lambda filename: suffix in filename and "muscle_spindles" in filename
muscle_spindles_files = ["out/{}".format(filename) for filename in sorted(list(np.extract(np.array(list(map(condition, filenames))), filenames)))]
muscle_spindles_files = muscle_spindles_files[::args.files_stride]


condition = lambda filename: suffix in filename and "motoneurons" in filename
motoneurons_files = ["out/{}".format(filename) for filename in sorted(list(np.extract(np.array(list(map(condition, filenames))), filenames)))]
motoneurons_files = motoneurons_files[::args.files_stride]

print("Number of input files: muscle spindles: {}, motoneurons: {}".\
  format(len(muscle_spindles_files), len(motoneurons_files)))


# load data
muscle_spindles_data = py_reader.load_data(muscle_spindles_files)
motoneurons_data = py_reader.load_data(motoneurons_files)

# create plots

if not args.plot_frequency_progression:
  fig,axes = plt.subplots(2,2,figsize=(12,6),sharex=True, gridspec_kw={'width_ratios': [10, 1]})
else:
  fig,axes = plt.subplots(3,2,figsize=(12,8),sharex=True, gridspec_kw={'width_ratios': [10, 1], "height_ratios": [1,1,0.5]})

# set global parameters for font sizes
plt.rcParams.update({
  'lines.linewidth': 2,
  'lines.markersize': 10,
  'lines.markeredgewidth': 2,
  'font.size': 10,
  'axes.labelsize': 10,
  'xtick.labelsize': 10,
  'ytick.labelsize': 10,
})

# ---------------------
# plot muscle spindles
component_name_input = "(P)modell/L"
component_name_output = "modell/primary_afferent"

t_values = None
values_output = None
values_input = None

# loop over datasets at different times
for i,dataset in enumerate(muscle_spindles_data):
  
  # get the data for the current timestep
  data_input = py_reader.get_values(dataset, "parameters", component_name_input)
  data_output = py_reader.get_values(dataset, "algebraics", component_name_output)
  
  if data_input is None:
    print("No data found for muscle spindles or component '{}' does not exist.\n".format(component_name_input))
  if data_output is None:
    print("No data found for muscle spindles or component '{}' does not exist.\n".format(component_name_output))
  
  # create arrays the first time
  if values_output is None:
    values_input = np.zeros((len(data_input), len(muscle_spindles_data)))   # each column is the data for one timestep, for multiple neurons
    values_output = np.zeros((len(data_output), len(muscle_spindles_data)))   # each column is the data for one timestep, for multiple neurons
    t_values = np.zeros((len(muscle_spindles_data)))
  
  # store values
  values_input[:,i] = data_input
  values_output[:,i] = data_output
  t_values[i] = dataset['currentTime']

# restrict data to specified end_time
start_index = None
for i in range(len(t_values)):
  if t_values[i] >= args.start_time and start_index is None:
    start_index = i
  if t_values[i] > args.end_time:
    end_index = i
    t_values = t_values[start_index:end_index]
    values_input = values_input[:,start_index:end_index]
    values_output = values_output[:,start_index:end_index]
    break
if args.end_time == np.inf:
  t_values = t_values[start_index:]
  values_input = values_input[:,start_index:]
  values_output = values_output[:,start_index:]
    
# plot lines for all timesteps
# loop over neurons
n = values_output.shape[0]
axes[0,0].plot([0,max(t_values)],[1,1], ":", color=(0.5,0.5,0.5))
for i in range(n):
  color = next(axes[0,0]._get_lines.prop_cycler)['color']
  axes[0,0].plot(t_values, values_input[i,:], ':', color=color)
  if i == 0:
    ax2 = axes[0,0].twinx()
  ax2.plot(t_values, values_output[i,:], '-', color=color)

# set title and axis labels
axes[0,0].set_title('Muscle spindles (number: {})'.format(n))
axes[0,0].set_ylabel('Sensed muscle stretch\n(dotted lines)', fontsize=12)
axes[0,0].grid(axis="x")
ax2.set_ylabel('Spindle response [mV]\n(solid lines)', fontsize=12)

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
    values_input = np.zeros((len(data_output), len(motoneurons_data)))   # each column is the data for one timestep, for multiple neurons
    values_output = np.zeros((len(data_output), len(motoneurons_data)))   # each column is the data for one timestep, for multiple neurons
    t_values = np.zeros((len(motoneurons_data)))
  
  # store values
  values_input[:,i] = data_input
  values_output[:,i] = data_output
  t_values[i] = dataset['currentTime']

# restrict data to specified end_time
start_index = None
for i in range(len(t_values)):
  if t_values[i] >= args.start_time and start_index is None:
    start_index = i
  if t_values[i] > args.end_time:
    end_index = i
    print("start_index: {}, end_index: {}".format(start_index, end_index))
    t_values = t_values[start_index:end_index]
    values_input = values_input[:,start_index:end_index]
    values_output = values_output[:,start_index:end_index]
    break
if args.end_time == np.inf:
  print("start_index: {}".format(start_index))
  t_values = t_values[start_index:]
  values_input = values_input[:,start_index:]
  values_output = values_output[:,start_index:]

if args.plot_firing_times:
  # ------
  # plot only firing times
  
  n_motor_units = values_output.shape[0]
  end_time = np.max(t_values)
  firing_threshold = np.mean(values_output)*3
  print("firing threshold: {}".format(firing_threshold))
  
  firing_times_mu = []
  for mu_no in range(n_motor_units):
    times = []
    current_value_active = False
    for i,value in enumerate(values_output[mu_no,:]):
      if value > firing_threshold and not current_value_active:
        current_value_active = True
        times.append(t_values[i])
      else:
        current_value_active = False
    firing_times_mu.append(times)
  #axes[1].plot(t_values, values_output[i,:], '-', color="k")
  
  print("end time: {}".format(end_time))
  
  # determine frequencies
  frequencies_mus = []
  for times in firing_times_mu:
    
    f = 1000 / (times[1]-times[0])
    frequencies = [f]
    for i in range(1,len(times)-1):
      f = 2000 / (times[i+1]-times[i-1])
      frequencies.append(f)
    frequencies.append(f)
    
    frequencies = np.array(frequencies)
    frequencies_mus.append(frequencies)
  
  min_frequency = np.min([np.min(m) for m in frequencies_mus])
  max_frequency = np.max([np.max(m) for m in frequencies_mus])
  
  norm = matplotlib.colors.Normalize(vmin=min_frequency,vmax=max_frequency)
  fig.colorbar(cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.viridis), ax=axes[1,1], pad=0.1, panchor=(0.0, 0.5))
      
  # plot line for each MU
  for (mu_no,times) in enumerate(firing_times_mu):
    
    axes[1,0].plot([0,end_time],[mu_no+1,mu_no+1],"k-")

    # plot stimulation times as vertical bar
    bar_height = 0.2
    
    color_values = (frequencies_mus[mu_no] - min_frequency) / (max_frequency-min_frequency)
    colors = matplotlib.cm.viridis(color_values)
    
    for i,t in enumerate(list(times)):
      axes[1,0].plot(t, mu_no+1, marker="+", color=colors[i])
    
    #axes[1].plot(list(times),[mu_no+1 for _ in times],marker="+",color=colors)
  
  # set title and axis labels
  axes[1,0].set_title('Motor neurons (number: {})'.format(n_motor_units))
  #axes[1,0].set_xlabel('time [ms]', fontsize=12)
  axes[1,0].set_ylabel('Motor units [-]', fontsize=12)
  axes[1,0].grid(axis="x")
  
  if args.plot_frequency_progression:
    l = []
    for mu_no in range(len(firing_times_mu)):
      for i in range(len(firing_times_mu[mu_no])):
        l.append((frequencies_mus[mu_no][i], firing_times_mu[mu_no][i]))
    freq_times = sorted(l, key=lambda a: a[1])
    
    fs = [x[0] for x in freq_times]
    ts = [x[1] for x in freq_times]
    
    ts = np.array([ts]).reshape(-1, 1)
    fs = np.array([fs]).reshape(-1, 1)
    
    #print("fs: {}".format(fs))
    #print("ts: {}".format(ts))
    
    # GPR, hyperparameters such as kernel width will be optimized
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, RBF, ConstantKernel
    
    kernel = ConstantKernel(constant_value=1.0) + RBF(length_scale=800, length_scale_bounds=(800.0, 1000.0))
    gpr = GaussianProcessRegressor(kernel=kernel).fit(ts, fs)
    
    frequency_values = gpr.predict(np.array([t_values]).reshape(-1,1))
    #print("frequency_values: {}".format(frequency_values))
    
    axes[2,0].plot(t_values, frequency_values, '-', label="With spindles")
    
    if args.load_frequency_file != "":
      with open(args.load_frequency_file, "rb") as f:
        (other_t_values,other_frequency_values) = pickle.load(f)
        
        n_values = min(len(other_frequency_values), len(t_values))
        
        axes[2,0].plot(other_t_values[0:n_values], other_frequency_values[0:n_values], '--', label="Without spindles")
        
    if args.store_frequency_file != "":
      with open(args.store_frequency_file, "wb") as f:
        print("Dump frequency values to file \"{}\"".format(args.store_frequency_file))
        pickle.dump((t_values,frequency_values), f)
    
    axes[2,0].set_ylim(10,max_frequency)
      
    # set title and axis labels
    axes[2,0].set_xlabel('Time [ms]', fontsize=12)
    axes[2,0].set_ylabel('Firing frequency [Hz]', fontsize=12)
    axes[2,0].legend()
    axes[2,0].grid(axis="x")
  
    axes[2,1].set_visible(False)
    
else:

  # ------
  # plot lines for all timesteps
  # loop over neurons
  n = values_output.shape[0]
  for i in range(n):
    color = next(axes[1,0]._get_lines.prop_cycler)['color']
    axes[1,0].plot(t_values, values_output[i,:], '-', color=color)
    if i == 0:
      ax2 = axes[1,0].twinx()
    ax2.plot(t_values, values_input[i,:], ':', color=color)
    
  # set title and axis labels
  axes[1,0].set_title('Motor neurons (number: {})'.format(n))
  axes[1,0].set_xlabel('time [ms]', fontsize=12)
  axes[1,0].set_ylabel('voltage [mV]\n(solid lines)', fontsize=12)
  ax2.set_ylabel('input current [uA]\n(dotted lines)', fontsize=12)

axes[0,1].set_visible(False)
axes[1,1].set_visible(False)

# show plot window
plt.tight_layout()
plt.savefig("plot.png")
plt.savefig("plot.pdf")
plt.show()
