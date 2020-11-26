#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize stimulation log files.
#

import sys
import numpy as np

import csv
import collections
import copy
import os
import time
import pickle
import py_reader    # reader utility for opendihu *.py files

filename = "stimulation.log"

if len(sys.argv) > 1:
  filename = sys.argv[1]
	
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

fiber_times = {}
fiber_mu_nos = {}

mu_times = {}
end_time = 0

# load data
with open(filename) as f:
  for line in f:
    if '#' in line:
      continue
    values = line.split(";")
    mu_no = int(values[0])
    fiber_no = int(values[1])
    times = [float(v)/1000.0 for v in values[2:-1]]
    
    fiber_mu_nos[fiber_no] = mu_no
    fiber_times[fiber_no] = times
    
    if not mu_no in mu_times:
      mu_times[mu_no] = set()
    
    for item in times:
      if item > end_time:
        end_time = item
      mu_times[mu_no].add(item)

print("end time: {}, n motor units: {}, n fibers: {}".format(end_time, len(mu_times), len(fiber_times)))


# plot firing times for mus
# --------------------------------
# prepare plot
fig = plt.figure()
# plot line for each MU
for (mu_no,times) in mu_times.items():
  
  plt.plot([0,end_time],[mu_no,mu_no],"k-")

  # plot stimulation times as vertical bar
  bar_height = 0.2
  plt.plot(list(times),[mu_no for _ in times],'k+')
  #for time in times:
  #  plt.plot([time,time],[mu_no-bar_height,mu_no+bar_height],'b-')

plt.xlabel('time [s]')
plt.ylabel('motor unit no [-]')

# plot firing times for fibers
# --------------------------------
fig = plt.figure()
# plot line for each fiber
for (fiber_no,times) in fiber_times.items():
  
  # plot stimulation times as vertical bar
  bar_height = 0.2
  plt.plot(list(times),[fiber_no for _ in times],'+')
  #for time in times:
  #  plt.plot([time,time],[mu_no-bar_height,mu_no+bar_height],'b-')

plt.xlabel('time [s]')
plt.ylabel('fiber no [-]')
plt.title('firing times for individual fibers units')

# plot firing times for fibers on mus
# --------------------------------
fig = plt.figure()

# plot line for each MU
for (mu_no,times) in mu_times.items():
  plt.plot([0,end_time],[mu_no,mu_no],"k-")

for (fiber_no,times) in fiber_times.items():
  
  mu_no = fiber_mu_nos[fiber_no]
  # plot stimulation times as vertical bar
  bar_height = 0.2
  plt.plot(list(times),[mu_no for _ in times],'+')
  #for time in times:
  #  plt.plot([time,time],[mu_no-bar_height,mu_no+bar_height],'b-')

plt.xlabel('time [s]')
plt.ylabel('motor unit no [-]')
plt.title('firing times for motor units')

plt.show()

#fig.patch.set_visible(False)
#ax.axis('off')



