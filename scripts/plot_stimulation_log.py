#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Script to visualize stimulation log files, possibly also with DEMUSE output and GRU prediction output.
# If used with GRU output, you have to manually adjust the parameters of the study in line ca. 288.
#
# usage:
# plot_stimulation_log.py [<stimulation.log> [<matlab file> [<gru file>]]]
#
# By default it looks for files "stimulation.log", "emg_results.mat" and "ipt_prediction.csv" if not specified as command line argument.
# The file "ipt_prediction.csv" has to contain the IPT as semicolon separated times, one line for every MU.
#

import sys
import numpy as np

import csv
import collections
import copy
import os
import time
import pickle
import traceback

filename = "stimulation.log"
matlab_filename = "emg_results.mat"
gru_filename = "ipt_prediction.csv"

# specify a plot time range, only works if no additional predicted data is loaded
begin_plot_time = 0
end_plot_time = np.inf

print("usage: plot_stimulation_log.py [<stimulation.log> [<matlab file> [<gru file>]]]\n")

# parse command line args
if len(sys.argv) > 1:
  filename = sys.argv[1]
if len(sys.argv) > 2:
  matlab_filename = sys.argv[2]
if len(sys.argv) > 3:
  gru_filename = sys.argv[3]
  
if begin_plot_time != 0 or end_plot_time != np.inf:
  matlab_filename = None
  gru_filename = None
  print("Note, specifying a plot range does not work with matlab or gru files!")
else:
  print("\nSearching for files \"{}\", \"{}\", \"{}\".".format(filename,matlab_filename,gru_filename))
  
# import needed packages from matplotlib
show_plot = True
if not show_plot:
  import matplotlib as mpl
  mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator

fiber_times = {}
fiber_mu_nos = {}

mu_times = {}
matlab_mus = None
end_time = 0

# load stimulation.log data
with open(filename) as f:
  for line in f:
    if '#' in line:
      continue
    values = line.split(";")
    
    # parse columns
    mu_no = int(values[0])
    fiber_no = int(values[1])
    times = [float(v)/1000.0 for v in values[2:-1] if begin_plot_time <= float(v)/1000.0 < end_plot_time]
    
    fiber_mu_nos[fiber_no] = mu_no
    fiber_times[fiber_no] = times
    
    # add empty set in mu_times for key mu_no
    if not mu_no in mu_times:
      mu_times[mu_no] = set()
    
    # add firing times in mu_times for key mu_no
    for item in times:
      if item > end_time:
        end_time = item
      mu_times[mu_no].add(item)

print("End time: {}, number of motor units: {}, number of fibers: {}".format(end_time, len(mu_times), len(fiber_times)))

def parse_mu_predictions(mu_firing_times, filename, find_mu_no=True, input_data_index_to_mu_no=None):
  """
  Parse given firing times of several MUs, compute timeshift compared to ground thruth and RoA metric
  :param mu_times: a list of lists with firing times
  :param filename: only for debugging output
  :return: a dict with key=MU no and value=dict with "times","onset","frequency","endtime","rate_of_agreement"
  """

  mu_frequencies_matlab = {}
  mu_onset_time_matlab = {}
  mu_endtime_matlab = {}
  mu_starttime_matlab = {}
  data_index_to_mu_no = {}
  
  # iterate over the found mus in the matlab file
  for i,firing_times in enumerate(mu_firing_times):
    
    # determine the firing frequency
    period_times = firing_times[1:] - firing_times[0:-1]
    frequency = 1/np.median(period_times)
    mu_frequencies_matlab[i] = frequency
    
    # determine the onset time, i.e., first real firing time, discard one outlier
    mu_onset_time_matlab[i] = np.min(firing_times)
    
    if (firing_times[1] - mu_onset_time_matlab[i]) > 4*frequency:
      print("outlier: {}, {}, {}, {}".format(firing_times[0], firing_times[1], firing_times[2], firing_times[3]))
      mu_onset_time_matlab[i] = firing_times[1]
    
    mu_endtime_matlab[i] = np.max(firing_times)
    mu_starttime_matlab[i] = np.min(firing_times)
    
  # map the MUs from the matlab file to the MUs from ground truth
  new_mu_times_matlab = {}
  new_mu_onset_time_matlab = {}
  new_mu_frequencies_matlab = {}
  matlab_mus = {}
  n_discarded_mus = 0
  
  # iterate over input data
  n_found_mus = len(mu_firing_times)
  for i in range(n_found_mus):
    
    # determine which mu_no the current mu no. i is
    if find_mu_no:
      scores = []
      mu_nos = []
      for mu_no in mu_times.keys():
        # score solely on onset time
        score1 = (mu_onset_time_ground_truth[mu_no] - mu_onset_time_matlab[i])**2   
        score2 = (mu_frequencies_ground_truth[mu_no] - mu_frequencies_matlab[i])**2
        score = score1 + score2
        scores.append(score)
        mu_nos.append(mu_no)

      if np.min(scores) > 9:
        print("Discard MU no. {} from {} file with onset time {}s".format(i, filename, mu_onset_time_matlab[i]))
        n_discarded_mus += 1
        continue

      mu_no = mu_nos[np.argmin(scores)]
      data_index_to_mu_no[i] = mu_no
    else:
      #mu_no = input_data_index_to_mu_no[i]   # for gCKC trained GRU
      mu_no = i+1                            # for stimulation.log trained GRU
    
    matlab_mus[mu_no] = {
      "original_data_index": i,
      "times": mu_firing_times[i],
      "onset": mu_onset_time_matlab[i],
      "frequency": mu_frequencies_matlab[i],
      "starttime": mu_starttime_matlab[i],
      "endtime": mu_endtime_matlab[i],
    }

  n_found_mus = len(matlab_mus)

  # compute timeshifts and metrics
  for mu_no in matlab_mus.keys():
    
    # compute timeshift
    distances = []
    for predicted_time in matlab_mus[mu_no]["times"]:
      ground_truth_times = list(mu_times[mu_no])
      index = np.argmin([abs(predicted_time - t) for t in ground_truth_times])
      distance_to_nearest_match = ground_truth_times[index] - predicted_time
      distances.append(distance_to_nearest_match)
    matlab_mus[mu_no]["timeshift"] = np.median(distances)

    tolerance = 5e-3
    tolerance = 10e-3

    # compute FP and TP metrics
    n_false_positives = 0
    for predicted_time in matlab_mus[mu_no]["times"]:
      ground_truth_times = list(mu_times[mu_no])
      index = np.argmin([abs(predicted_time+matlab_mus[mu_no]["timeshift"] - t) for t in ground_truth_times])
      distance = predicted_time+matlab_mus[mu_no]["timeshift"] - ground_truth_times[index]
      if abs(distance) > tolerance:
        #print("  offset: {} -> FP".format(distance))
        n_false_positives += 1
      #else:
      #  print("  offset: {} -> TP".format(distance))
    matlab_mus[mu_no]["n_false_positives"] = n_false_positives
    matlab_mus[mu_no]["n_true_positives"] = len(matlab_mus[mu_no]["times"]) - n_false_positives

    # compute FN and RoA metric
    n_false_negatives = 0
    for firing_time in mu_times[mu_no]:
      if firing_time < matlab_mus[mu_no]["starttime"]:
        continue
      if firing_time > matlab_mus[mu_no]["endtime"]:
        break
      predicted_times = matlab_mus[mu_no]["times"]
      index = np.argmin([abs(firing_time - (t+matlab_mus[mu_no]["timeshift"])) for t in predicted_times])
      distance = predicted_times[index]+matlab_mus[mu_no]["timeshift"] - firing_time
      if abs(distance) > tolerance:
        n_false_negatives += 1
    matlab_mus[mu_no]["n_false_negatives"] = n_false_negatives
    matlab_mus[mu_no]["precision"] = (float)(matlab_mus[mu_no]["n_true_positives"]) / (matlab_mus[mu_no]["n_true_positives"] + matlab_mus[mu_no]["n_false_positives"])
    matlab_mus[mu_no]["rate_of_agreement"] = (float)(matlab_mus[mu_no]["n_true_positives"]) / (matlab_mus[mu_no]["n_true_positives"] + matlab_mus[mu_no]["n_false_positives"] + matlab_mus[mu_no]["n_false_negatives"])
    print(" MU {} from {}: FP: {}, TP: {}, FN: {} -> RoA: {}".format(mu_no, filename, matlab_mus[mu_no]["n_false_positives"], matlab_mus[mu_no]["n_true_positives"], matlab_mus[mu_no]["n_false_negatives"], matlab_mus[mu_no]["rate_of_agreement"]))

  return matlab_mus, n_found_mus, n_discarded_mus, data_index_to_mu_no

# parse ground truth
mu_frequencies_ground_truth = {}
mu_onset_time_ground_truth = {}
  
# iterate over ground truth mu times
for (mu_no,times) in mu_times.items():
  times = list(sorted(times))
  if len(times) < 2:
    break
  
  # determine the firing frequency
  period_times = np.array(times[1:]) - np.array(times[0:-1])
  frequency = 1/np.median(period_times)
  mu_frequencies_ground_truth[mu_no] = frequency
  
  # determine the onset time, i.e., first real firing time, discard one outlier
  mu_onset_time_ground_truth[mu_no] = np.min(times)
  
  if (times[1] - mu_onset_time_ground_truth[mu_no]) > 4*frequency:
    mu_onset_time_ground_truth[mu_no] = times[1]
    
matlab_mus = {}
gru_mus = {}

# load matlab data
if matlab_filename is not None:
  import scipy.io
  try:
    print("\nTry to load and process data from matlab file \"{}\"".format(matlab_filename))
    matlab_data = scipy.io.loadmat(matlab_filename)
    # matlab_data["MUPulses"][n_found_mus][i]   # firing times 
    # matlab_data["IPTs"][n_found_mus][t]       # signal how probable it is that a MU fires at the time t, threshold with 0.1 to obtain labels
    
    
    # iterate over the found mus in the matlab file
    mu_times_matlab = []
    n_found_mus = matlab_data["MUPulses"].shape[1]
    for i in range(n_found_mus):
      
      # store the firing times
      mu_times_matlab.append(np.array([t/2000.0 for t in matlab_data["MUPulses"][0][i][0]]))
      
    # determine time shifts, map to MU no. and compute rate of agreement metric
    matlab_mus, n_found_mus_matlab, n_discarded_mus_matlab, data_index_to_mu_no \
      = parse_mu_predictions(mu_times_matlab, "matlab")
    
  except:
    print("Could not load data from a matlab file \"{}\".".format(matlab_filename))
    traceback.print_exc()
    pass
    
# load GRU predicted data
if gru_filename is not None:
  try:
    print("\nTry to load and process data from GRU file \"{}\"".format(gru_filename))
    
    # load mu firing times from gru filename
    mu_times_gru = []
    with open(gru_filename, "r") as f:
      
      # load line of file
      lines = f.readlines()
      for line in lines:
        
        # convert to values
        if ";" in line:
          firing_times = line.split(";")
          firing_times = map(float, firing_times)
        else:
          continue
        
        # transform input data from ms to s and add offset for test data
        firing_times = np.array([t/2000.0 for t in firing_times])
        
        # predicted data starts at data point 60400 with sampling frequency of 2kHz, i.e. at time 30.2s
        firing_times = np.array([t + 30+0.2 for t in firing_times])
        print("load IPT in time range range [{},{}]s".format(float(np.min(firing_times)), float(np.max(firing_times))))
        
        
        #firing_times = np.array([t for t in firing_times])
        mu_times_gru.append(firing_times)

    #print("mu_times_gru: {}".format(mu_times_gru))

    # determine time shifts, map to MU no. and compute rate of agreement metric
    gru_mus, n_found_mus_gru, n_discarded_mus_gru, _ = parse_mu_predictions(mu_times_gru, "gru", False, data_index_to_mu_no)
    
  except:
    print("Could not load data from a GRU file \"{}\".".format(gru_filename))
    traceback.print_exc()
    pass
    
# debugging output
for (mu_no,times) in mu_times.items():
  
  print("MU {:2} ground truth, onset: {:.3f} s, frequency: {:.2f} Hz".format(mu_no, mu_onset_time_ground_truth[mu_no], mu_frequencies_ground_truth[mu_no]))
  if mu_no in matlab_mus.keys():
    d = matlab_mus[mu_no]
    print("             matlab onset: {:.3f} s, frequency: {:.2f} Hz (timeshift: {:.5f}s), RoA: {:.3f}, data end: {} s".
      format(d["onset"], d["frequency"], d["timeshift"], d["rate_of_agreement"],  d["endtime"]))
  
  if mu_no in gru_mus.keys():
    d = gru_mus[mu_no]
    print("                gru onset: {:.3f} s, frequency: {:.2f} Hz (timeshift: {:.5f}s), RoA: {:.3f}, precision: {:.3f}, data end: {} s".
      format(d["onset"], d["frequency"], d["timeshift"], d["rate_of_agreement"], d["precision"],  d["endtime"]))
    

if matlab_mus:
  print("Matlab file \"{}\" contains {} matching MUs, {} discarded.".format(matlab_filename, n_found_mus_matlab, n_discarded_mus_matlab))

if gru_mus:
  print("GRU file \"{}\" contains {} matching MUs, {} discarded.".format(gru_filename, n_found_mus_gru, n_discarded_mus_gru))

# set global parameters for font sizes
plt.rcParams.update({'font.size': 18})
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 10

# plot firing times for mus
# --------------------------------
# prepare plot
#fig = plt.figure(figsize=(6,3))
fig = plt.figure(figsize=(12,6))
# plot line for each MU
for (mu_no,times) in mu_times.items():
  
  plt.plot([0,end_time],[mu_no,mu_no],"k-")

  # plot stimulation times as vertical bar
  bar_height = 0.2
  plt.plot(list(times),[mu_no for _ in times],'k+')
  #for time in times:
  #  plt.plot([time,time],[mu_no-bar_height,mu_no+bar_height],'b-')

# plot matlab predictions
if matlab_mus:
  for (mu_no,d) in matlab_mus.items():
    plt.plot(list(d["times"]+d["timeshift"]),[mu_no+0.2 for _ in d["times"]],'r+')
    #plt.plot(list(d["times"]),[mu_no+0.4 for _ in d["times"]],'r+')
  
# plot gru predictions
if gru_mus:
  for (mu_no,d) in gru_mus.items():
    plt.plot(list(d["times"]+d["timeshift"]),[mu_no-0.2 for _ in d["times"]],'b+')
  
ax = fig.gca()
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
#ax.set_xlim(15.9,18.5)
#ax.set_ylim(13.7,18.7)
#plt.subplots_adjust(bottom=0.22,left=0.2)
plt.xlabel('Time [s]')
plt.ylabel('Motor unit no [-]')

# plot firing times for fibers
# --------------------------------
if end_plot_time - begin_plot_time < 10:
  fig = plt.figure()
  # plot line for each fiber
  for (fiber_no,times) in fiber_times.items():
    
    # plot stimulation times as vertical bar
    bar_height = 0.2
    plt.plot(list(times),[fiber_no for _ in times],'+')
    #for time in times:
    #  plt.plot([time,time],[mu_no-bar_height,mu_no+bar_height],'b-')

  ax = fig.gca()
  ax.yaxis.set_major_locator(MaxNLocator(integer=True))
  plt.xlabel('Time [s]')
  plt.ylabel('Fiber no [-]')
  plt.title('firing times for individual fibers')

plt.show()

#fig.patch.set_visible(False)
#ax.axis('off')



